#!/usr/bin/env python3
desc="""Get alternative poly-A sites from BAM files and report their frequencies and counts.
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 28/2/2022
"""

import csv, gzip, os, sys
import matplotlib.pyplot as plt, numpy as np, pysam, seaborn as sns 
from datetime import datetime
from scipy.signal import find_peaks
from pybedtools import BedTool

def get_peaks(counts, prominence=5, dist=10, min_count=25):
    """Return peaks and their boundaries"""
    peaks, props = find_peaks(counts, prominence=prominence)
    if len(peaks)<1: return [], []
    centers, proms = [peaks[0]], [props['prominences'][0]]
    se  = [[props['left_bases'][0], props['right_bases'][0]],]
    for c, p, s, e in zip(peaks[1:], props['prominences'][1:], 
                          props['left_bases'][1:], props['right_bases'][1:]):
        if s-dist<=se[-1][-1]:
            if e>se[-1][-1]: se[-1][-1] = e
            if p>proms[-1]: centers[-1], proms[-1] = c, p
        else:
            centers.append(c)
            se.append([s, e])
    # filter by minimum number of reads
    _centers, _se = [], []
    for c, (ps, pe) in zip(centers, se):
        if counts[ps:pe].sum() < min_count: continue
        _centers.append(c)
        _se.append((ps, pe))
    centers, se = np.array(_centers), np.array(_se)
    return centers, se

def get_counts(sams, ref, s, e, strand, mapq=10, firststrand=False, extend=0):
    """Return counts of read ends"""
    counts = np.zeros((len(sams), e-s+extend), dtype="int")
    for fi, sam in enumerate(sams): 
        for a in sam.fetch(ref, s, e): 
            if a.mapq<mapq or a.is_secondary or a.is_supplementary: continue
            # get transcript strand (cDNA is antisense, DRS is sense)
            is_reverse = not a.is_reverse if firststrand else a.is_reverse
            if strand=="+" and is_reverse or strand=="-" and not is_reverse: continue
            # for reverse algs store the beginning of the alignment
            p = a.pos-s+extend if is_reverse else a.aend-s
            # skip if end outside of gene body
            if p<0 or p>=counts.shape[1]: continue
            counts[fi, p] += 1
    return counts

def get_pA_sites(outfn, gff, fnames, samples=[], min_count=25, mapq=10,
                 firststrand=False, extend=2000, min_frac=0.05, 
                 verbose=False, logger=sys.stderr.write):
    """Report frequencies and counts for alternative polyA sites from BAM"""
    # get sample names
    if not samples: samples = [os.path.basename(fn)[:-4] for fn in fnames]
    else: assert len(samples)==len(fnames), "Number of samples has to match number of BAM files!"
    
    # open SAM files
    sams = [pysam.AlignmentFile(fn) for fn in fnames]
    # get gene locations
    ann = BedTool(gff)
    genes = ann.filter(lambda x: x[2].endswith("gene"))

    # prepare tsv file
    with gzip.open(outfn, 'wt') if outfn.endswith('gz') else open(outfn, 'wt') as out:
        tsv = csv.writer(out, delimiter='\t')
        header = ["chrom", "start", "end", "gene_name", "alt_site_no", "strand", "feature", "peak_center"]
        header += ["%s freq"%s for s in samples] + ["%s counts"%s for s in samples]
        tsv.writerow(header)
        j = k = 0
        for i, g in enumerate(genes, 1):
            # g.attrs 'logic_name': 'araport11'
            feature, name, ref, s, e, strand = g[2], g.name, g.chrom, g.start, g.end, g.strand
            if name.startswith('gene'): name = name[5:]
            if not i%10: logger(f"{i:,} {ref} {name} so far: {k:,} pA sites from {j:,} genes    \r")
            # get read ends
            counts = get_counts(sams, ref, s, e, strand, mapq, firststrand, extend)
            if counts.sum() < min_count: continue
            # get alternative pA sites
            centers, se = get_peaks(counts.sum(axis=0), min_count=min_count)
            if len(centers)<1: continue
            j += 1
            # make sure we always count from the last pA site
            if not firststrand and strand=="+" or firststrand and strand=="-":
                centers, se = centers[::-1], se[::-1]
                # we need to extend antisense transcripts
                s -= extend
            # get counts for each alt pA site
            peak_counts = np.array([counts[:, ps:pe].sum(axis=1).flatten() for ps, pe in se])
            sel = peak_counts.sum(axis=1)/peak_counts.sum(axis=1).max() >= min_frac
            # store frequencies and counts for each pA
            for ii, (center, (ps, pe)) in enumerate(zip(centers[sel], se[sel]), 1):
                c = counts[:, ps:pe].sum(axis=1).flatten()
                f = np.round(c/peak_counts.sum(axis=0), 3)
                tsv.writerow((ref, s+ps, s+pe, name, ii, strand, feature, s+center, *f, *c))
            k += ii
            #if i>100: break
    logger("Results for %s alternative pA sites from %s genes reported to %s\n"%(k, j, outfn))

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0a')   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-a", "--gff", required=1,
                        help="gene annotation in gtf/gff format")
    parser.add_argument("-b", "--bam", nargs="+",
                        help="input BAM files")
    parser.add_argument("-o", "--out", default="pA_sites.tsv.gz",
                        help="output file name [%(default)s]")
    parser.add_argument("-d", "--minDepth", default=25, type=int,
                        help="min depth of coverage (for each strand separately) [%(default)s]")
    parser.add_argument("-q", "--mapq", default=10, type=int,
                        help="min mapping quality [%(default)s]")
    parser.add_argument("-s", "--samples", nargs="*",
                        help="sample names [infer from BAM names]")
    parser.add_argument("--firststrand", action='store_true', 
                        help="antisense strand sequencing ie cDNA [sense strand sequencing ie DRS]")
    parser.add_argument("-e", "--extend", default=1000,  
                        help="extend gene by [%(default)s bases]")
    parser.add_argument("-f", "--min_frac", default=0.05,  
                        help="report pA sites with at least [%(default)s] reads of max pA site for given gene")
    
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    get_pA_sites(o.out, o.gff, o.bam, o.samples, o.minDepth, o.mapq,
                 o.firststrand, o.extend, o.min_frac, o.verbose)
  
if __name__=='__main__': 
    t0 = datetime.now()
    os.setpgrp() # create new process group, become its leader    
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except Exception as err:
        import signal, traceback
        sys.stderr.write(traceback.format_exc()+"\n")
        os.killpg(0, signal.SIGTERM) # terminate all processes in my group
    finally:
        dt = datetime.now()-t0
        sys.stderr.write("#Time elapsed: %s    \n"%dt)
