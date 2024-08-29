# polyTailor

- [Brief description](#Brief-description)
- [Installation](#Installation)
- [How to run?](#How-to-run?)
- [Examples](#Examples)
- [Test dataset](#Test-dataset)
- [Citation](#Citation) 
- [Issues](#Issues)

## Brief description 
PolyTailor is a software to study RNA tails.
It predicts: 
* per-read polyA tail lengths estimates
* per-read tail heterogeneity (i.e. non-A features)
 
This software is meant to be used with Nano3P-seq cDNA libraries,
and can work with Nano3P-seq libraries sequenced using R9 and R10 flowcells.

## Installation

You'll need:
- [samtools](https://github.com/samtools/htslib/) v1.19+
- [minimap2](https://github.com/lh3/minimap2/) v2.28+
- [IsoQuant](https://github.com/ablab/IsoQuant) 3.5+
- Python 3.8+ with following packaged installed: matplotlib numpy parasail pybedtools pysam pandas scipy seaborn htseq

All above can be installed using conda and pip: 

```bash
conda create -c conda-forge -c bioconda -n polyTailor python=3.10 samtools minimap2 isoquant
conda activate polyTailor
pip install matplotlib numpy parasail pybedtools pysam pandas scipy seaborn htseq
mkdir -p ~/src && cd ~/src
git clone https://github.com/novoalab/polyTailor.git
```

- [dorado](https://github.com/nanoporetech/dorado) v0.7.2+ and basecalling models. 
Here, we assume Linux x64 system. For other systems, please see
[dorado page](https://github.com/nanoporetech/dorado?tab=readme-ov-file#installation)
.

```bash
mkdir -p ~/src && cd ~/src
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.7.2-linux-x64.tar.gz
tar xpfz dorado-0.7.2-linux-x64.tar.gz
echo 'export PATH=~/src/dorado-0.7.2-linux-x64/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

dorado download --directory ~/src/dorado/models --model dna_r9.4.1_e8_hac@v3.3
dorado download --directory ~/src/dorado/models --model dna_r9.4.1_e8_sup@v3.3
dorado download --directory ~/src/dorado/models --model dna_r10.4.1_e8.2_400bps_hac@v5.0.0
dorado download --directory ~/src/dorado/models --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0
```


## How to run?

PolyTailor can be executed as follows (steps 3-4 are optional): 

### 1. Basecall (and demultiplex) your reads saving the `mv` table in BAM file
```bash
dorado basecaller -x cuda:all --emit-moves -r MODEL [--kit-name BARCODING_KIT] pod5_dir > reads.bam
```

For the most accurate poly-T composition calling we recommend using the latest `sup` 
[model](https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models).
If barcoding `--kit-name` is provided, barcode will be reported in `barcode` column. 


### 2. Align the reads to the genome passing `dorado` tags to the resulting BAM file

```bash
samtools fastq -F2304 -T mv,ts,BC reads.bam|minimap2 -y -ax splice:hq genome.fa -|samtools sort --write-index -o algs.bam
```

### 3. Annotate alternative transcript ends

```bash
src/get_transcript_ends.py --firststrand -q0 -o transcript_ends.tsv.gz -a genome.gtf -b algs.bam [algs2.bam ... algsN.bam]
```

### 4. Associate reads to transcripts

```bash
isoquant.py --complete_genedb --data_type nanopore -o isoquant -r genome.fa -g genome.gtf --stranded reverse --bam algs.bam
```

### 5. Estimate poly-T tail length and composition

```bash
src/get_pT.py -o pT.tsv.gz -b algs.bam \
  -e transcript_ends.flt.tsv.gz.bed
  -i <(zgrep -v '^#' isoquant/OUT/OUT.read_assignments.tsv.gz | cut -f1,4,6,9)
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
   - the pT is immediately following primer (otherwise `not_continuous`)
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

## Test dataset

You can find test data and example outputs in [test](/test). 

## Citation

If you find this work useful, please cite: 
Begik O*, Pryszcz LP*, Niazi AM, Valen E, Novoa EM.
Nano3P-seq: a protocol to chart the coding and non-coding transcriptome at single molecule resolution.
(in preparation)

## Issues

If you have an issue running this code, please open a new Github issue.
Please take a look at previous issues, even if closed. Thanks! 
