#!/usr/bin/env python3
desc="""Report poly-T tail length from basecalled (BAM with `ts` and `mv` tags) N3PS library. 

Primer sequence (-p) should include terminal -TTT
in order to facilitate better primer alignment. 

TODO:
- composition score
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Torredembarra/Barcelona/Mizerów, 9/08/2024
"""

import os
import sys
import pysam
import numpy as np
import parasail
import pysam
import re
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
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

def get_pt_boundaries(seq, pt_start, pat=re.compile(r'((T{2,}\w?)+)', re.IGNORECASE), before=10, after=50):
    start = end = 0
    for m in pat.finditer(seq[pt_start-before:pt_start+after]):
        s, e = m.span()
        if s>2*before: continue
        s = pt_start + s - before + 1 # we include nonT start/end
        e = pt_start + e - before - 1 # we include nonT start/end
        if e-s > end - start:
            start, end = s, e
    return start-1, end

def alg2pt(a, primer, profile, open_penalty, extend_penalty, extension=5, analyse_first=250):
    """Return sequence flanking pt, pt start, pt end, pt length, pt sequence, alignment score, identity and cigar"""
    seq = a.get_forward_sequence()
    r = parasail.sw_trace_striped_profile_16(profile, seq[:analyse_first], open_penalty, extend_penalty) # 41.6 µs vs 13.4 µs :200
    # filter TODO: scoring
    #if r.score<50: continue
    cigar = r.cigar.decode.decode() # 2.62 µs
    score, ts, te, qs, qe = r.score, r.cigar.beg_ref, r.end_ref, r.cigar.beg_query, r.end_query
    cigar, ts, qs, m, nm, qalg, talg = correct_alg(cigar, ts, qs) # 2.6 µs
    # get polyT length in bases
    end_skip = len(primer)-qe-1
    if qs: cigar = f"{qs}S" + cigar
    if end_skip: cigar += f"{end_skip}S"
    pt_start = te + end_skip - 2
    # get pt boundaries
    s, e = get_pt_boundaries(seq, pt_start)
    # unload tags
    tags = {k: v for k, v in a.tags}
    identity = round(1-nm/qalg, 3)
    mv = np.array(tags['mv'])
    stride, move = mv[0], mv[1:]
    move_pos = np.argwhere(move==1).flatten()
    steps = move_pos[1:]-move_pos[:-1]
    per_base = steps[e+extension+1:].mean()
    pt_len = round(steps[s-extension:e+extension+1].sum()/per_base - 2*extension, 1)
    #return seq[pt_start-10:pt_start+50], s, e, pt_len, seq[s:e], r.score, identity, cigar
    return pt_len, per_base, s, e, seq[s-10:s], seq[s:e], r.score, identity, cigar

def get_pT(out, bam, readidsfn, primer, #mapq=10, max_dist=10, 
           minscore=20, scoring=(2, 3, 3, 2), logger=sys.stderr):

    # init primer aligner
    match_score, mismatch_penalty, open_penalty, extend_penalty = scoring
    matrix = parasail.matrix_create("ACGT", match_score, -mismatch_penalty)
    profile = parasail.profile_create_16(primer, matrix)

    # load readids
    readids = {}
    if readidsfn:
        if logger: logger.write(f"Loading read ids...\n")
        for l in open(readidsfn, "rt"):
            if not l[:-1]: continue
            ldata = l[:-1].split("\t")
            readids[ldata[0]] = "\t".join(ldata[1:])
        if logger: logger.write(f" {len(readids):,} entries loaded.\n")
    
    if logger: logger.write(f"Processing reads...\n")
    sam = pysam.AlignmentFile(bam)
    #ref2len = {r: l for r, l in zip(sam.references, sam.lengths)
    comment = ""
    out.write("read_id\tpt_length\tper_base\tpt_start\tpt_end\tbefore_pt\tpt_seq\tscore\tidentity\tcigar\tcomments\n")
    for ai, a in enumerate(sam, 1):
        if not ai%1000: logger.write(f" {ai:,} \r")
        readid = a.qname
        if readids:
            if readid not in readids: continue
            comment = readids[readid]
        #e = a.aend if a.is_reverse else a.pos
        if a.is_secondary: continue # or a.is_supplementary or a.mapq<mapq or abs(e-ref2len[ref])>max_dist: continue
        pt_data = alg2pt(a, primer, profile, open_penalty, extend_penalty)        
        out.write(f"%s\t%s\t%.1f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(readid, *pt_data, comment))
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
    parser.add_argument("-i", "--readids", default="", 
                        help="file with read ids [use all reads]")
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("w"), 
                        help="output name [stdout]")
    parser.add_argument("-p", "--primer",
                        default="CAGCACCT"+"ACTTGCCTGTCGCTCTATCTGCAGAGCAGAG"+"TTT",
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

    get_pT(o.out, o.bam, o.readids, o.primer)
    
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
