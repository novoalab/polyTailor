#!/usr/bin/env python3
desc="""Report poly-T tail length from basecalled (BAM with `ts` and `mv` tags) N3PS library. 

Primer sequence (-p) should include terminal -TTT
in order to facilitate better primer alignment. 

TODO:
- composition score
- report single pT estimate per read - using terminal alg
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Torredembarra/Barcelona/Mizerów, 9/08/2024
"""

import gzip
import os
import sys
import pysam
import numpy as np
import parasail
import pysam
import re
import traceback
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import HTSeq
from datetime import datetime

def correct_alg(cigar, ts, qs, w=6, k=3, pat=re.compile(r'\d+')):
    """Correct alignment skipping initial deletion
    
    Return corrected: 
    - cigar 
    - target start
    - query start
    - number of anchors
    - edit distance (NM)
    - aligned query bases
    """
    idx = n_anchors = nm = qalg = talg = 0
    pe = 0
    diff = w-k
    for i, m in enumerate(pat.finditer(cigar)):
        e = m.end()
        l = int(cigar[pe:e])
        op = cigar[e:e+1]
        pe = e+1
        # skip initial deletion
        if i==0 and op in "ID":
            idx = pe
            if op == "D": ts += l
            elif op == "I": qs += l
        # update number of anchors
        elif op=="=":
            n = (l-diff) // k
            if n>0:
                n_anchors += n
            qalg += l
            talg += l
        # count mismatches
        elif op!="S":
            nm += l
            if op!="D": qalg += l
            if op!="I": talg += l
    return cigar[idx:], ts, qs, n_anchors, nm, qalg, talg

def get_pt_boundaries(seq, pt_start, max_start=10, before=0, 
                      pat=re.compile(r'((T{2,}\w?)+)', re.IGNORECASE)):
    start = end = 0
    for m in pat.finditer(seq[pt_start-before:]):
        s, e = m.span()
        if s>max_start: continue
        if e-s > end - start:
            start = s + pt_start + before
            end = e + pt_start + before - 1 # last base is notT
    return start, end

def alg2pt(a, primer, profile, open_penalty, extend_penalty,
           extension=5, analyse_first=250, minscore=41):
    """Return sequence flanking pt, pt start, pt end, pt length, pt sequence, alignment score, identity and cigar"""
    info = "OK"
    empty = (0, 0, 0, "", "")
    # get size of clipped region of the read (from the 5' read end)
    clipped = 0
    if a.cigar:
        idx = -1 if a.is_reverse else 0
        clipped = a.cigar[idx][1] if a.cigar[idx][0]==4 else 0
    if clipped<len(primer):
        return "no_clip", *empty
    # get forward sequence of a read trimming aligned part
    seq = a.get_forward_sequence()[:min(analyse_first, clipped+2*extension)]
    # find primer in the read sequence
    r = parasail.sw_trace_striped_profile_16(profile, seq, open_penalty, extend_penalty) # 41.6 µs vs 13.4 µs :200
    # filter TODO: scoring
    #if r.score<50: continue
    cigar = r.cigar.decode.decode() # 2.62 µs
    score, ts, te, qs, qe = r.score, r.cigar.beg_ref, r.end_ref, r.cigar.beg_query, r.end_query
    if score<minscore:
        return "no_primer", *empty        
    cigar, ts, qs, m, nm, qalg, talg = correct_alg(cigar, ts, qs) # 2.6 µs
    #identity = round(1-nm/qalg, 3)
    # get polyT length in bases
    end_skip = len(primer)-qe-1
    if qs: cigar = f"{qs}S" + cigar
    if end_skip:
        cigar += f"{end_skip}S"
        info = "not_complete"
    pt_start = te + end_skip + 1 # - 2
    # get the longest pT stretch
    s, e = get_pt_boundaries(seq, pt_start-1, 2*extension)
    if e-s<2:
        return "no_pT", *empty
    if s>pt_start:
        info = "not_continuous"
    # unload tags
    tags = {k: v for k, v in a.tags}
    mv = np.array(tags['mv'])
    stride, move = mv[0], mv[1:]
    move_pos = np.argwhere(move==1).flatten()
    steps = move_pos[1:]-move_pos[:-1]
    # let's take per base from aligned bases?
    per_base = steps[:s-extension].mean() #steps[e+extension+1:].mean()
    # here we should capture and remove stalls
    pt_len = round(steps[s-extension:e+extension+1].sum()/per_base - 2*extension, 1)
    return info, pt_len, per_base, s, seq[s-10:s], seq[s:e]

def get_pT(out, bam, readidsfn, endsfn, primer, scoring=(2, 3, 3, 2), logger=sys.stderr):
    """
    """
    # create gz handle
    if out.name.endswith(".gz"):
        out.close()
        out = gzip.open(out.name, "wt")
    # init primer aligner
    match_score, mismatch_penalty, open_penalty, extend_penalty = scoring
    matrix = parasail.matrix_create("ACGT", match_score, -mismatch_penalty)
    profile = parasail.profile_create_16(primer, matrix)

    # load readids
    readids = {}
    ends = {}
    if readidsfn:
        if logger: logger.write(f"Loading read ids...\n")
        f = gzip.open(readidsfn, "rt") if readidsfn.endswith('.gz') else open(readidsfn, "rt")
        for l in f:
            if not l[:-1]: continue
            ldata = l[:-1].split("\t")
            readids[ldata[0]] = "\t".join(ldata[1:])
        f.close()
        if logger: logger.write(f" {len(readids):,} entries loaded.\n")
    if endsfn:
        ends = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        for f in HTSeq.BED_Reader(endsfn):
            ends[f.iv] = (f.name, f.thick.start)
    
    if logger: logger.write(f"Processing reads...\n")
    sam = pysam.AlignmentFile(bam, check_sq=False)
    #ref2len = {r: l for r, l in zip(sam.references, sam.lengths)
    out.write("read_id\tbarcode\tmapq\tfilter\tpt_length\tper_base\tpt_start\tbefore_pt\tpt_seq\ttranscript_end\tdistance\tcomment\n")
    outline = "%s\t%s\t%s\t%s\t%s\t%.1f\t%s\t%s\t%s\t%s\t%s\t%s\n"
    #try:
    for ai, a in enumerate(sam, 1):
        if not ai%1000: logger.write(f" {ai:,} \r")
        # skip not primary algs (2304) and unmapped reads (4)
        if a.flag&2308: continue
        readid = a.qname
        comment = alt_end = dist = ""
        # get barcode
        tags = {k: v for k, v in a.tags}
        barcode = tags["BC"] if "BC" in tags else "unknown"
        # store alternative transcript end info
        if ends and a.reference_name in ends.chrom_vectors:
            # --firststrand
            e, strand = (a.aend, "+") if a.is_reverse else (a.pos, "-")
            p = HTSeq.GenomicInterval(a.reference_name, e, e+1, strand)
            alt_end = ";".join([v[0] for iv, v in ends[p].steps() if v])
            dist = ";".join(map(str, [abs(e-v[1]) for iv, v in ends[p].steps() if v]))
        # store additional info ie isoquant
        if readids:
            if readid in readids:
                comment = readids[readid]
        pt_data = alg2pt(a, primer, profile, open_penalty, extend_penalty)
        out.write(outline%(readid, barcode, a.mapq, *pt_data, alt_end, dist, comment))
    #except Exception as e: sys.stderr.write("\n"+traceback.format_exc()+"\n")
    if logger: logger.write(f" {ai:,} reads processed.\n")

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0a') 
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-b", "--bam", required=True, 
                        help="input BAM file with `ts` and `mv` tags")
    parser.add_argument("-e", "--ends", default="", 
                        help="file with alternative ends")
    parser.add_argument("-i", "--readids", default="", 
                        help="file with read ids [use all reads]")
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("w"), 
                        help="output name [stdout]")
    parser.add_argument("-p", "--primer",
                        default="CAGCACCT"+"ACTTGCCTGTCGCTCTATCTGCAGAGCAGAG",
                        help="N3PS primer sequence [%(default)s]")
    '''
    parser.add_argument("-m", "--mapq", default=10, type=int, 
                        help="min. mapping quality [%(default)s]")
    parser.add_argument("--firststrand", action='store_true', 
                        help="antisense strand sequencing ie cDNA [sense strand sequencing ie DRS]")
    parser.add_argument("-c", "--complete", action='store_true', 
                        help="count only reads covering entire tRNA reference (+5 bases of each oligo)")
    '''
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    get_pT(o.out, o.bam, o.readids, o.ends, o.primer)
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
