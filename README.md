# PolyTailor

- [Brief description](#Brief-description)
- [Installation](#Installation)
- [How to run?](#How-to-run?)
- [Examples](#Examples)
- [Citation](#Citation) 
- [Issues](#Issues)

## Brief description 
PolyTailor is a software to study RNA tails.
It predicts: 
* per-read polyA tail lengths estimates
* per-read tail heterogeneity (i.e. non-A features)
 
This software is meant to be used with Nano3P-seq cDNA libraries, and can work with Nano3P-seq libraries sequenced using R9 and R10 flowcells.

## Installation

You'll need Python 3.8+ and install following packaged 
using [pip](https://pip.pypa.io/):

```bash
pip install matplotlib numpy parasail pybedtools pysam pandas scipy seaborn htseq
```

## How to run?

PolyTailor can be executed as follows (steps 3-4 are optional): 

### 1. Basecall (and demultiplex) your reads saving the `mv` table in BAM file 
using [dorado](https://github.com/nanoporetech/dorado).
```bash
dorado basecaller -x cuda:all --emit-moves -r MODEL [--kit-name BARCODING_KIT] pod5_dir > reads.bam
```
For the most accurate poly-T composition calling we recommend using the latest `sup` 
[model](https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models).
If barcoding `--kit-name` is provided, barcode will be reported in `barcode` column. 


### 2. Align the reads to the genome passing `dorado` tags to the resulting BAM file

```bash
samtools fastq -F2304 -T mv,ns,pt,ts,BC reads.bam|minimap2 -y -ax splice:hq genome.fa -|samtools sort --write-index -o algs.bam
```

### 3. Annotate alternative transcript ends

```bash
src/get_transcript_ends.py --firststrand -q0 -o transcript_ends.tsv.gz -a genome.gtf -b algs.bam [algs2.bam ... algsN.bam]
```

### 4. Associate reads to transcripts
using [IsoQuant](https://github.com/ablab/IsoQuant)

```bash
isoquant.py --complete_genedb --data_type nanopore -o isoquant -r genome.fa -g genome.gtf --stranded reverse --bam algs.bam
```

### 5. Estimate poly-T tail length and composition combining all above info

```bash
src/get_pT.py -o pT.tsv -b algs.bam [-e transcript_ends.tsv.gz -i <(zgrep -v '^#' isoquant/OUT/OUT.read_assignments.tsv.gz | cut -f1,4,6,9)]
```

This will produce a TAB-delimited file with following columns:
1. read_id
2. barcode - detected barcode (`unknown` if no barcode detected)
3. mapq - mapping quality
4. filter - `OK` means that following filters were passed
   - read has 5' clipped part corresponding to: adapter, N3PS primer and pT (otherwise `no_clip`)
   - N3PS primer was detected in the 5' clipped part (otherwise `no_primer`)
   - the N3PS primer aligned end-to-end (otherwise `not_complete`)
   - pT sequence was detected between primer and aligned transcript (otherwise `no_pT`)
   - the pT is immediately following primer (otherwise not_continuous)
5. pt_len - estimated poly-T length. 
6. per_base - helicase speed estimated from mv table (mean number of chunks per base)
7. pt_start - poly-T start in read sequence
8. before_pt - sequence before detected poly-T (terminal bases of poly-A)
9. pt_seq - poly-T sequence composition
10. transcript_end - 
11. distance - distance from the transcript end
12. comments - additional fields passed from `-i / --readids` file

Note, there may be multiple comments columns,
depending on provided `-i / --readids` file.

For example, for `isoquant` example above, you'll see:
12. isoform_id
13. assignment_type
14. additional_info


## Examples

### a) Old N3Pseq oligo `CAGCACCT CTTCCGATCACTTGCCTGTCGCTCTATCTTC` (used in original Nano3P-seq protocol, Nat Methods 2022)

- 30 bases long poly-A tails
```bash
read_id	barcode	mapq	filter	pt_length	per_base	pt_start	before_pt	pt_seq
66726005-d641-4c06-8b00-4621f94ec03f	unknown	60	OK	27.3	1.5	102	CTCTATCTTC	TTTTTTTTTTTTTTTTTTT
89456d91-eb3a-46ad-be21-63814a7cfe9b	unknown	9	OK	33.6	3.2	100	CTCTATCTTC	TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
ddfd08ad-c316-4197-acbb-6f5a00468397	unknown	60	OK	25.8	2.2	104	CTCTATCTTC	TTTTTTTTTTTTTTTTTTTTTTTTT
c7b689ae-9dd4-482b-b2b8-1e4edbc46a55	unknown	60	OK	23.3	2.2	109	CTCTATCTTC	TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
3e11c499-2a54-4373-a360-65fe0a2d0a9b	unknown	60	OK	37.7	1.8	100	CTCTATCTTC	TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
```

- 30 bases long poly-A tails with internal G bases
```bash
read_id	barcode	mapq	filter	pt_length	per_base	pt_start	before_pt	pt_seq
00dbf8ad-5efe-41ec-817c-d67ea271e454    unknown 51      OK      12.8    5.0     106     CTCTATCTTC      TTTTCTTTTCTTTTCTTTTTTTTTTTT
6f44ebc0-5128-40c6-b72d-599367f6ab1f    unknown 60      OK      21.5    2.3     104     CTCTATCTTC      TTTTCTTTTCTTTTCTTTTTTTTTT
fc395607-5f52-415d-92f8-62bf6f22628a    unknown 54      OK      22.0    2.4     109     CTCTATCTTC      TTTTCTTTTCTTTTCTTTTTTT
468f9a11-4941-4644-bfa2-dd788c3184ff    unknown 3       OK      25.5    2.2     108     CTCTATCTTC      TTTTCTTTTCTTTTCTTTTTTT
b1c908b0-0ad8-41fd-9e65-83182e09af76    unknown 60      OK      31.6    2.2     100     CTCTATCTTC      TTTTCTTTTCTTTTCTTTTTTTTTTTT
```

- 30 bases long poly-A tails with a single terminal U base `-U`
```bash
read_id	barcode	mapq	filter	pt_length	per_base	pt_start	before_pt	pt_seq
9c2dc8f5-2520-4216-bb6a-7b63a9811a97	unknown	60	not_continuous	26.6	2.2	105	TCTATCTTCA	TTTTTT
2c10382e-8cca-4be5-842c-582e6c90ec1d	unknown	60	not_continuous	20.3	2.0	104	TATATCTTCA	TTTTTTT
1c1d96c5-d64d-48c8-88b2-a9efa10859ff	unknown	60	not_continuous	32.3	1.8	105	TCTATCTTCA	TTTTTTTTTTTTTTTTTTTTTTT
197b468e-20e5-4b0b-9526-c72e32911256	unknown	60	not_continuous	33.3	2.1	100	TCTATCTTCA	TTTTTTTTTTTTTTTT
798987a3-5b67-4dba-b737-be9ff84b753b	unknown	60	not_continuous	25.8	2.3	104	TCTATCTTCA	TTTTTTTTTTT
```

- 30 bases long poly-A tails with 3 terminal U bases `-UUU`
```bash
read_id	barcode	mapq	filter	pt_length	per_base	pt_start	before_pt	pt_seq
89894967-58ff-415d-b007-44119d12eada	unknown	3	not_continuous	10.5	3.3	107	ATCTTAGAAA	TTTTTTTTTTTTTTTTTTTTT
67430d4b-7ec1-4a4f-a92f-bd103d39d32d	unknown	60	not_continuous	18.9	2.4	108	TATCTTCAAA	TTTTT
ec8cf495-5c2c-4b87-b507-a887c533a8e8	unknown	34	not_continuous	14.6	2.9	105	TATCTTCAAA	TTTTTTTTTTTTTTTTTTTT
848fb138-282c-4baa-8932-f325d6778184	unknown	60	not_continuous	24.9	1.9	105	TATCTTCAAA	TTTTTTTTTTTTTTTTTTTTTTTTTTTTT
d277c17f-d690-494f-b858-bed07c959c89	unknown	60	not_continuous	17.6	2.6	111	TATCTTCAAA	TTTTTTTTTTTTTTTTTTTTTTT
```

- 30 bases long poly-A tails with 5 terminal U bases `-UUUUU`
```bash
read_id	barcode	mapq	filter	pt_length	per_base	pt_start	before_pt	pt_seq
20b1a13a-6eab-4d0b-9d10-0d71c58a7d38	unknown	60	not_continuous	35.9	1.9	111	TCTTCAAAAA	TTTTTTTTTTTTTTTTTTTTTTTTTT
1e4c7b8f-be74-4f1c-8dc1-50c0c6cd4c60	unknown	60	not_continuous	22.1	2.1	109	ATCTTCAAAA	TTTTTTTTTTTTTT
fcd0dfe8-3fac-4e62-ade8-c064b46e5f87	unknown	60	not_continuous	16.9	2.3	111	TCTTCAAAAA	TTTTTTTTTTTTTTTTT
bec9bb54-aa7a-4492-ba9c-00f0f3f729fe	unknown	60	not_continuous	24.1	2.5	113	TCTTCAAAAA	TTTTTTTTTTTTTTT
8e2e89d7-1b79-46e1-a585-e02d2ab30a2e	unknown	60	not_continuous	12.3	3.6	106	TCTTCAAAAA	TTTTTTTTTTTTTTTTTTTTTT
```

- 30 bases long poly-A tails with 5 terminal C bases `-CCCCC`
```bash
read_id	barcode	mapq	filter	pt_length	per_base	pt_start	before_pt	pt_seq
f0f7fe2b-3f7c-4474-bad9-0721557a35d3	unknown	60	not_continuous	22.0	1.9	110	TCTTCGGGGG	TTTTTTTTTT
2ffc2b3c-59d7-45db-b994-6c3140a1570c	unknown	54	not_continuous	24.0	1.7	111	TGTTCGGGGG	TTTTTTTTTTTTTTTTTTTTTT
904389b4-3830-4147-ac78-1a5f030b6690	unknown	60	not_continuous	19.1	2.2	108	TCTTCGGGGG	TTTTTTTTTTTTTTTTTTT
4196161c-959d-432a-bfca-d9d14f3acc90	unknown	60	not_continuous	15.6	2.5	110	TCTTCGGGGG	TTTTTTTTTTTTTTTTTTTTT
7d0d7548-a714-4efa-bc02-958db3740d8e	unknown	60	not_continuous	18.1	2.4	108	TCTTCGGGGG	TTTTTTTTTTTTTTTTTTTTT
```

- 30 bases long poly-A tails with 5 terminal G bases `-GGGGG`
```bash
read_id	barcode	mapq	filter	pt_length	per_base	pt_start	before_pt	pt_seq
d64295f5-4e5d-4b5c-b24d-29b2606318b4	unknown	35	not_continuous	21.1	1.6	122	TGCCCCCCCC	TTTTTTT
af7fa273-9fd3-4ad9-a1f9-9691b49e19be	unknown	60	not_continuous	22.9	2.1	110	TCTTCCCCCC	TTTTTTT
073bab93-2b78-4c7d-b1ff-34f3b67a394f	unknown	60	not_continuous	5.6	3.7	105	TCTTCCCCCC	TTTTTTTT
a1eea540-b878-4bdd-9251-9ac65f475200	unknown	60	not_continuous	19.0	2.6	107	CTTCCCCCCC	TTTTTTT
692be019-144c-486a-b329-81c0a44f93d8	unknown	60	not_continuous	27.5	2.5	101	ATCTTCCCCC	TTTTTTTTT
```

### b)  New N3Pseq oligo `CAGCACCT ACTTGCCTGTCGCTCTATCTGCAGAGCAGAG` (this work)

- yeast total RNA
```bash
read_id	barcode	mapq	filter	pt_length	per_base	pt_start	before_pt	pt_seq	transcript_end	distance	comment
7d7576b3-3872-4c22-baf0-620f35f0917f	unknown	7	OK	57.8	2.1	105	CCAGAGCAAG	TTTT
11cdb20c-1808-4a5d-a6c6-7850daafb4b8	unknown	60	OK	36.3	2.5	100	CAGAGCAGAG	TTTTTTTTTTTTTTTTTT	GDH3|1|protein_coding	17
dc61cb45-7fa1-4c09-af85-21ce91e28879	unknown	60	OK	35.6	3.5	106	CAGAGCAGAG	TTTTTTTTTTTTTTTTTTT	GDH3|1|protein_coding	18
19625470-110e-4315-8354-79fd7574aadd	unknown	0	no_pT	0	0.0	0	RDN25-1|2|rRNA	1	
5517395a-bc83-471c-8f03-867b46333fcc	unknown	0	no_pT	0	0.0	0	RDN25-1|2|rRNA	1	
```

## Citation

If you find this work useful, please cite: 
Begik O*, Pryszcz LP*, Niazi AM, Valen E, Novoa EM. Nano3P-seq: a protocol to chart the coding and non-coding transcriptome at single molecule resolution. (in preparation)

## Issues

If you have an issue running this code, please open a new Github issue. Please take a look at previous issues, even if closed. Thanks! 
