# Test

We provide following examples:
- [precomputed](#precomputed): aligned BAM files with `mv` table
  (and optionally precomputed transcript ends and read-transcript table)
- [raw reads](#raw-reads): pod5 files, reference genome FastA and GTF

## Precomputed

We provide precomputed datasets for: 
- DNA standards - synthetic cDNA oligos complementary to RNA strand terminated with
  - `*A`: A-only tail with varying lengths (0, 15, 30, 60, 90 or 120);  
  - `*U`, `5G`, `5C`: 30 bases long poly-A tail that ends by either Us (1, 3 or 5), Gs (5) or Cs (5); 
  - `IntG`: 15 As followed by 3 GAAAA repeats
- Yeast total RNA sequenced with R9 (old N3PS oligo) and R10 (new N3PS oligo) chemistry

Above data is already basecalled and aligned,
so you'll perform only the last step of polyTailor. 


### DNA standards

Note, DNA standards were sequenced using old N3Pseq oligo `CAGCACCT CTTCCGATCACTTGCCTGTCGCTCTATCTTC`. 

We prefiltered only the reads for which
the alignment starts in expected reference posittion (+/-1).
For example, for 0A, 60A and 120A,
the alignments should start at position
1, 61 and 121 of the reference, respectively 
(adding 12 bases that are typically lost from 5'-end by ONT sequencing). 

First download the data

```bash
cd ~/src/polyTailor/test
wget https://public-docs.crg.es/enovoa/public/lpryszcz/src/polyTailor/test/{ref,minimap2} -q --show-progress -r -c -nc -np -nH --cut-dirs=6 --reject="index.html*"
```

Then process all aligned reads with polyTailor:
```bash
cd test
for f in minimap2/DNA_standards/*.bam; do
  echo $f;
  python ../src/get_pt.py -b $f -o $f.pT.tsv.gz \
    -p CAGCACCTCTTCCGATCACTTGCCTGTCGCTCTATCTTC;
done

for f in minimap2/DNA_standards/*.bam.pT.tsv.gz; do echo $f; zcat $f|head; done
```

This will produce an output file for each BAM file similar to this:

- minimap2/DNA_standards/30A.bam.pT.tsv.gz: 30 bases long poly-A tails

```bash
read_id	barcode	mapq	filter	pt_length	per_base	primer_end	pt_start	before_pt	pt_seq	transcript_end	distance	comment
a5befef8-9780-4d0a-88f8-27941bb13c5d	unknown	60	OK	26.3	2.2	108	108	CTCTATCTTC	TTTTTTTTTTTTTTTTTTTTTTTTTTTT
3b41842a-fbb4-4c84-bde7-e3c1adb645c3	unknown	60	OK	22.7	2.4	106	106	CTCTATCTTC	TTTTTTTTTTTTTTTTTTTTTTTTTTT
8632b9b8-0c4d-4f3c-82e3-8142704869c4	unknown	60	OK	22.3	1.9	80	80	CTCTATCTTC	TTTTTTTTTTTTTTTTTTTTTTTT
7a5748d3-c905-48fe-8b16-8c886c7c5bfc	unknown	60	OK	39.5	2.0	103	103	CTCTATCTTC	TTTTTTT
b5f88ac9-924d-4d64-bdba-234c6ac03a7d	unknown	60	not_complete	23.1	2.2	100	99	ACTCTGTCTT	TTTTTTTTTTTTTTTTTTTTTTTT
b016d2c7-bf3c-4352-bc9c-a1ceeb05dc4b	unknown	60	no_primer	0	0.0	0	0
25056f68-1672-4476-9aa3-da3d36571150	unknown	60	not_continuous	33.4	2.0	105	106	TCTATTACTC	TTT
5fa7b265-720f-4baa-9c5f-c0aff976e89f	unknown	60	not_continuous	26.7	2.0	111	112	TCTATCTTTC	TTTTTTTTTTTTTTTTTT
80a1a47e-a6f8-41dd-810c-ffc4e8382ed7	unknown	60	OK	18.9	2.3	99	99	CTCTATCTTC	TTTTTTTTTTTTTTTTTTTTTTT
```

Additional figures can be generated using
[/notebook/DNA_standards.ipynb](/notebook/DNA_standards.ipynb)

![DNA standards: length](/notebook/DNA_standards.length.png "DNA standards: length")
![DNA standards: composition](/notebook/DNA_standards.composition.png "DNA standards: composition")

### Yeast total RNA

Here, new N3PS oligo `CAGCACCT ACTTGCCTGTCGCTCTATCTGCAGAGCAGAG` was used,
therefore we don't need to change `-p`.


```bash

```

- yeast total RNA
```bash
read_id	barcode	mapq	filter	pt_length	per_base	pt_start	before_pt	pt_seq	transcript_end	distance	comment
7d7576b3-3872-4c22-baf0-620f35f0917f	unknown	7	OK	57.8	2.1	105	CCAGAGCAAG	TTTT
11cdb20c-1808-4a5d-a6c6-7850daafb4b8	unknown	60	OK	36.3	2.5	100	CAGAGCAGAG	TTTTTTTTTTTTTTTTTT	GDH3|1|protein_coding	17
dc61cb45-7fa1-4c09-af85-21ce91e28879	unknown	60	OK	35.6	3.5	106	CAGAGCAGAG	TTTTTTTTTTTTTTTTTTT	GDH3|1|protein_coding	18
19625470-110e-4315-8354-79fd7574aadd	unknown	0	no_pT	0	0.0	0	RDN25-1|2|rRNA	1	
5517395a-bc83-471c-8f03-867b46333fcc	unknown	0	no_pT	0	0.0	0	RDN25-1|2|rRNA	1	
```


## Raw reads

We provide raw datasets for: 
- Yeast total RNA sequenced with R9 (old N3PS oligo) and R10 (new N3PS oligo) chemistry

### Yeast total RNA

If you wish to repeat the full polyTailor analysis descibed in the manuscript,
follow these steps: 

0. Download reference genome and raw reads in pod5 format

```bash
cd ~/src/polyTailor/test
wget https://public-docs.crg.es/enovoa/public/lpryszcz/src/polyTailor/test/{ref,reads} -q --show-progress -r -c -nc -np -nH --cut-dirs=6 --reject="index.html*"
```

If you have previous results in `minimap2/yeast`,
either remove or rename this directory.

```
mv minimap2/yeast minimap2/yeast.old
```


1. Basecall


- R9 run (0.5h for `hac` vs 1.5h for `sup`)

```bash
cd ~/src/polyTailor/test
mkdir -p dorado/yeast
dorado basecaller \
  ~/src/dorado/models/dna_r9.4.1_e8_hac@v3.3 \
  reads/yeast/R9 -r \
  --kit-name EXP-NBD104 \
  --no-trim \
  --emit-moves  \
  -x cuda:all > dorado/yeast/R9.bam
```

- R10 run (0.5h for `hac` vs 12h for `sup`)

```bash
dorado basecaller \ 
  ~/src/dorado/models/dna_r10.4.1_e8.2_400bps_hac@v5.0.0 \
  reads/yeast/R10 -r \
  --kit-name SQK-NBD114-24 \
  --no-trim \
  --emit-moves  \
  -x cuda:all > dorado/yeast/R10.bam
```

2. Align

```bash
mkdir -p minimap2/yeast
for f in dorado/yeast/*.bam; do
  echo `date` $f; 
  samtools fastq -T mv,ts,BC $f \
  | minimap2 -y -ax splice:hq -G2k ref/yeast.fa.gz - \
  | samtools sort --write-index -o minimap2/yeast/$(basename $f);
done
```

3. Transcript ends

```bash
../src/get_transcript_ends.py --firststrand -q 0 -e 1000 \
  -a ref/yeast.gtf.gz -b minimap2/yeast/*.bam \
  -o minimap2/yeast/transcript_ends.tsv.gz 
```

4. Assign reads to transcripts

```bash
isoquant.py --complete_genedb --data_type nanopore --stranded reverse \
  -r ref/yeast.fa.gz -g ref/yeast.gtf.gz --bam minimap2/yeast/*.bam \
  -o isoquant
zgrep -v '^#' isoquant/OUT/OUT.read_assignments.tsv.gz \
  | cut -f1,4,6,9 | gzip > minimap2/yeast/read_assignments.tsv.gz
```

5. Estimate transcript ends and composition

- R9 (old N3PS oligo)

```bash
f=minimap2/yeast/R9.bam; 
../src/get_pT.py -o $f.pT.tsv.gz -b $f \
  -e minimap2/yeast/transcript_ends.flt.tsv.gz.bed \
  -i minimap2/yeast/read_assignments.tsv.gz \
  -p CAGCACCTCTTCCGATCACTTGCCTGTCGCTCTATCTTC
```

- R10 (new N3PS oligo)

```bash
f=minimap2/yeast/R10.bam; 
../src/get_pT.py -o $f.pT.tsv.gz -b $f \
  -e minimap2/yeast/transcript_ends.flt.tsv.gz.bed \
  -i minimap2/yeast/read_assignments.tsv.gz
```




