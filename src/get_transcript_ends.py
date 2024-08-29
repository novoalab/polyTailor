#!/usr/bin/env python3
desc="""Get alternative transcript ends (poly-A sites) from BAM files
and report their frequencies and counts.
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona/Mizerów, 9/08/2024
"""

import csv, gzip, os, sys, traceback
import matplotlib.pyplot as plt, numpy as np, pysam, seaborn as sns 
from datetime import datetime
from pybedtools import BedTool
from multiprocessing import Pool
from scipy.signal import find_peaks

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
        try:
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
        except Exception as e:
            sys.stderr.write("\n"+traceback.format_exc()+"\n")
    return counts

def _init_args(*args):
    global sams, mapq, firststrand, extend
    fnames, mapq, firststrand, extend = args
    sams = [pysam.AlignmentFile(fn) for fn in fnames]
    # ignore invalid
    np.seterr(invalid='ignore') #divide, over, under, invalid

def worker(args):
    feature, name, ref, s, e, strand = args
    counts = get_counts(sams, ref, s, e, strand, mapq, firststrand, extend)
    return counts

def get_pA_sites(outfn, gtf, fnames, samples=[], min_count=25, mapq=10,
                 firststrand=False, extend=2000, min_frac=0.05, threads=6, 
                 verbose=False, logger=sys.stderr.write):
    """Report frequencies and counts for alternative polyA sites from BAM"""
    if os.path.isfile(outfn):
        logger(f"File exists: {outfn}!\n")
        return
    # get sample names
    if not samples: samples = [os.path.basename(fn)[:-4] for fn in fnames]
    else: assert len(samples)==len(fnames), "Number of samples has to match number of BAM files!"
    
    logger(f"Loading genes...\n")
    # get gene locations
    ann = BedTool(gtf)
    genes = ann.filter(lambda x: x[2].endswith("gene"))
    # fill args
    args, refs = [], set()
    for i, g in enumerate(genes, 1):
        # g.attrs 'logic_name': 'araport11'
        ref, s, e, strand = g.chrom, g.start, g.end, g.strand
        feature = None
        for f in ("gene_biotype", "gene_type"):
            if f in g.attrs:
                feature = g.attrs[f]
                break
        if not feature: continue
        name = g.name # g.attrs 'gene_id' 'gene_name'==g.name
        args.append((feature, name, ref, s, e, strand))
        refs.add(ref)
    # add chrs without any genes as additional genes (ie sequin, CC etc)
    with pysam.AlignmentFile(fnames[0]) as sam:
        args += [("ext", r, r, 0, l, "+")
                 for r, l in zip(sam.references, sam.lengths) if r not in refs]
    logger(f"Processing {len(args):,} genes...\n")
    # open SAM files
    p = Pool(min(threads, len(args)), initializer=_init_args,
             initargs=(fnames, mapq, firststrand, extend))

    # prepare tsv file
    with gzip.open(outfn, 'wt') if outfn.endswith('gz') else open(outfn, 'wt') as out:
        tsv = csv.writer(out, delimiter='\t')
        header = ["chrom", "start", "end", "name", "usage", "strand", "peak_center"]
        header += ["%s freq"%s for s in samples] + ["%s counts"%s for s in samples]
        tsv.writerow(header)
        j = k = 0
        for i, counts in enumerate(p.imap(worker, args), 1):
            feature, name, ref, s, e, strand = args[i-1]
            if not i%10: logger(f"{i:,} {ref} {name} so far: {k:,} pA sites from {j:,} genes    \r")
            # get read ends
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
            peak_counts = np.array([counts[:, ps:pe].sum(axis=1).flatten() for ps, pe in se]); print(centers, se, peak_counts)
            sel = peak_counts.sum(axis=1)/peak_counts.sum(axis=1).max() >= min_frac
            # store frequencies and counts for each pA
            for ii, (center, (ps, pe)) in enumerate(zip(centers[sel], se[sel]), 1):
                c = counts[:, ps:pe].sum(axis=1).flatten()
                f = np.round(c/peak_counts.sum(axis=0), 3)
                usage = np.round(np.nanmean(c/peak_counts.sum(axis=0)), 3)
                _name = f"{name}|{ii}|{feature}"
                tsv.writerow((ref, s+ps, s+pe, _name, usage, strand, s+center, *f, *c))
            k += ii
            #if i>100: break
    logger("Results for %s alternative pA sites from %s genes reported to %s\n"%(k, j, outfn))
    # close pool
    p.close()
    p.join()

def filter_predictions(fn, logger=sys.stderr.write):
    """Keep only one prediction for pA sites shared by multiple genes by:
    - skip intronic (not implemented)
    - and assign to the gene with fewer pA sites (or shorter if the same number of pA sites)

    TBD:
    - work with overlaps (currently only start position
    - report shorter if two genes have the same amount of pA sites
    """
    def report_lines(out, outbed, pos2data, gene2count, gene2len={}):
        if not pos2data: return
        for lines in pos2data.values():
            # sort by gene with lower gene2count
            lines = sorted(lines, key=lambda x: gene2count[x[0]])
            # report only 1 line for overlap
            line = lines[0][1]
            out.write(line)
            # write bed
            if line.startswith('chrom\tstart'): continue
            ldata = line[:-1].split('\t')
            usage, strand, center = ldata[4:7]
            usage, center = float(usage), int(center)
            rgb = ",".join(map(str, np.round(usage*np.array((255 if strand=="-" else 0, 255 if strand=="+" else 0, 0)), 0).astype('int')))
            outbed.write("\t".join(ldata[:6]) + f"\t{center}\t{center+1}\t{rgb}\n")

    if fn.endswith(".gz"):
        out = gzip.open(fn[:-7]+".flt.tsv.gz", "wt")
        handle = gzip.open(fn, "rt")
    else:
        out = open(fn[:-4]+".flt.tsv", "wt")
        handle = open(fn, "rt")
    outbed = open(out.name+".bed", "wt")
    if logger: logger(f"Saving filtered predictions to: {out.name}\n")
    pchrom = ""
    pos2data = {} # this can be improved using array with pointer to line list
    gene2count = {}
    for l in handle:
        ldata = l[:-1].split('\t')
        chrom, s, e, name, score, strand = ldata[:6]
        gene = name.split("|")[0]
        if chrom!=pchrom:
            report_lines(out, outbed, pos2data, gene2count)
            pchrom = chrom
            pos2data = {}
            gene2count = {}
        k = f"{s}{strand}"
        if k in pos2data: pos2data[k].append((gene, l))
        else: pos2data[k] = [(gene, l)]
        if gene not in gene2count: gene2count[gene] = 1
        else: gene2count[gene] += 1
    report_lines(out, outbed, pos2data, gene2count)
    out.close(); outbed.close(); handle.close()

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0a')   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-a", "--gtf", required=1,
                        help="gene annotation in gtf format")
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
    parser.add_argument("-e", "--extend", default=1000, type=int,  
                        help="extend gene by [%(default)s bases]")
    parser.add_argument("-f", "--min_frac", default=0.05, type=float,  
                        help="report pA sites with at least [%(default)s] reads of max pA site for given gene")
    parser.add_argument("-t", "--threads", default=6, type=int, 
                        help="no. of threads to use [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    # get pA sites for all genes
    get_pA_sites(o.out, o.gtf, o.bam, o.samples, o.minDepth, o.mapq,
                 o.firststrand, o.extend, o.min_frac, o.threads, o.verbose)
    
    # pA sites with multi-genes: skip intronic and assign to the gene with fewer pA sites
    filter_predictions(o.out)
  
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
