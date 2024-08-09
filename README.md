# N3PS_pt

## Installation

You'll need to install several Python packages using pip:

```bash
pip install matplotlib numpy parasail pysam pandas seaborn
```

## How to run?
You can estimate poly-T tail length and composition using:
```bash
src/get_pT.py -b bam_with_move_table.bam > pT.tsv
```

Since, N3Pseq captures all transcripts (also fragmented ones),
you may want to subset the analysis only to reads that end at
one of the annotated poly-A sites.

DETAILS FROM OZ HERE!

Once you have a BED file, you can extract read_id, gene_name and feature_type using:
```bash
cut -f4,16-17 complete_reads.bed > read.ids
```
and use it with `get_pT.py` as follows:
```bash
src/get_pT.py -b bam_with_move_table.bam -i read.ids > pT.flt.tsv
```


This will produce a TAB-delimited file with following columns:
1. read_id
2. pt_len - estimated poly-T length. 
3. per_base - helicase speed estimated from mv table (mean number of chunks per base)
4. pt_start - poly-T start in read sequence
5. pt_end - poly-T end in read sequence
6. before_pt - sequence before detected poly-T (terminal bases of poly-A)
7. pt_seq - poly-T sequence composition
8. score - N3Pseq primer alignment score
9. identity - N3Pseq primer identity
10. cigar - N3Pseq primer alignment cigar
11. comments - additional fields passed from `-i / --readids` file

Note, there may be multiple comments columns,
depending on provided `-i / --readids` file. 


### Example outputs

1. Old N3Pseq oligo `CAGCACCT CTTCCGATCACTTGCCTGTCGCTCTATCTTC TTT`

- 30 bases long poly-A tails
```bash
read_id pt_length       per_base        pt_start        pt_end  before_pt       pt_seq  score   identity        cigar   comments
ebd3774d-8d97-4934-8b68-450a07ff93e5    34.3    2.2     105     115     TCGCTCTATC      TTCTTTTTTT      74      0.952   11=1I1=1X28=
95d75aca-4dbc-47cb-adda-77bc780c6b16    28.1    2.4     101     111     TCGCTCTATC      TTCTTTTTTT      84      1.0     42=
2f559f2b-51fc-436b-b9c9-cddb996a6a03    34.4    1.9     94      100     TCGCTCTATC      TTCTTT  79      0.976   19=1I22=
fc30b14f-e814-4ea5-bcc1-c3d0e40868be    34.5    1.8     100     130     TCGCTCTATC      TTCTTTTTTTTTTTTTTTTTTTTTTTTTTT  84      1.0     42=
5104d2f0-9ef3-401b-a9ef-7ba2bc37cb8c    35.7    1.6     112     118     GCACTCTATC      TTCTTT  30      0.727   9S3=1D1=1X2=1X1=1D3=1I3=1X1=2D1=1X13=
```

- 30 bases long poly-A tails with internal G bases
```bash
read_id pt_length       per_base        pt_start        pt_end  before_pt       pt_seq  score   identity        cigar   comments
c561064b-8ee3-44c9-8454-6b5810dd1f78    15.3    2.8     98      132     TCGCTCTATC      TTCTTTTTTCTTTTCTTTTCTTTTTTTTTTTTTT      84      1.0     42=
c8dac82d-84ab-4e20-8009-2167b7e3e663    28.3    1.8     100     126     TCGCTCTATC      TTCTTTTCTTTTCTTTTCTTTTTTTT      84      1.0     42=
e467dfa0-b598-483f-a8ae-cfd2d9004c0f    28.9    2.1     98      126     TCGCTCTATC      TTCTTTTCTTTTCTTTCTTTTTTTTTTT    84      1.0     42=
447d09fa-b8d5-475b-ba19-229ce1f2e766    24.9    2.4     91      119     TCGCTCTATC      TTCTTTTCTTTTCTTTTCTTTTTTTTTT    84      1.0     42=
3a8e11e7-238c-434f-bfb3-a37e6858ab4b    44.8    1.3     104     132     TCGCTCTATC      TTCTTTTCTTTTCTTTTCTTTTTTTTTT    84      1.0     42=
```

- 30 bases long poly-A tails with a single terminal U base (-U)
```bash
read_id pt_length       per_base        pt_start        pt_end  before_pt       pt_seq  score   identity        cigar   comments
4fc31a31-aaf9-481b-a9e2-cb8f04622485    35.0    2.4     108     129     CTATCTTTAA      TTTTTTTTTTTTTTTTTTTTT   71      0.974   19=1I18=4S
ebbac215-516c-480a-8bb4-833912240e2e    50.7    2.2     94      103     TCTATCTTCA      TTTTTTGTT       81      0.976   39=1D3=
13732215-ffbb-40e1-89a3-65e422e0d183    22.5    2.2     103     127     TCTATCTTCA      TTTTTTTTTTTTTTTTTTTTTTTT        81      0.976   39=1D3=
e41cc3d2-8855-4d0d-85a2-849fe9d193ff    17.3    2.7     104     123     GTATTTTTCA      TTTTTTTTTTTTTTTTTTT     46      0.786   22=2I1=1D2=2X2=1X2=1D1=1X3=1D3=
1142c841-8716-4956-90e1-372a3fc1a675    22.9    2.1     107     116     TCTATCTTCA      TTTTTTTTT       68      0.905   8=1D3=2X26=1D3=
```

- 30 bases long poly-A tails with 3 terminal U bases (-UUU)
```bash
read_id pt_length       per_base        pt_start        pt_end  before_pt       pt_seq  score   identity        cigar   comments
1ba6eb7c-2319-4c79-aa05-c967a13eeeb1    25.0    1.9     106     126     TATCTTCAAA      TTTTTTTTTTTTTTTTTTTT    78      1.0     39=3S
952c5aec-da57-41c9-9ab9-fd2c0c0f8e2d    29.3    1.7     108     133     TCTATTTAAA      TTTTTTTTTTTTTTTTTTTTTTTTT       71      0.974   35=1I2=4S
6fce380f-d07f-43c9-bbd4-4388b2ebbefb    20.6    1.8     106     127     ATCTCTCAAA      TTTTTTTTTTTTTTTTTTTTT   66      0.925   29=2I5=1I3=2S
917cf7cc-fcbc-41d7-a382-1a4e336b1940    22.1    2.1     106     127     TATCTTCAAA      TTTTTTTTTTTTTTTTTTTTT   78      1.0     39=3S
2d2d3344-00f8-4245-9015-cbcc80eb3d2d    25.8    1.9     115     139     ATCTTAGAAA      TTTTTTTTTTTTTTTTTTTTTTTT        76      1.0     38=4S
```

- 30 bases long poly-A tails with 5 terminal U bases (-UUUUU)
```bash
read_id pt_length       per_base        pt_start        pt_end  before_pt       pt_seq  score   identity        cigar   comments
0c6a6d34-d8fd-4f73-8f89-230d94e03992    31.6    1.9     104     124     TCTTCAAAAA      TTTTTTTTTTTTTTTTTTTT    57      0.872   9=1I17=1X1=1D1=1D1=1X7=3S
b532612e-8e89-4d3b-a0c3-7841d585bb1b    33.2    2.0     107     130     TCTTCAAAAA      TTTTTTTTTTTTTTTTTTTTTTT 75      0.974   8=1D31=3S
a8d4c6d5-91d8-4bd2-b875-93d0ea42415e    32.0    2.3     109     120     TCTTCAAAAA      TTTTTTTTTTT     78      1.0     39=3S
8e7f20dd-7b6b-4b6c-8178-3963557a60fd    23.5    2.1     119     133     TCTTCAAAAA      TTTTTTTTTTTTTT  78      1.0     39=3S
f62e13f2-0371-4c1e-9b46-3c9e04cebc9c    22.1    1.9     102     114     TCTTCAAAAA      TTTTTTTTTTTT    78      1.0     39=3S
```

- 30 bases long poly-A tails with 5 terminal C bases (-CCCCC)
```bash
read_id pt_length       per_base        pt_start        pt_end  before_pt       pt_seq  score   identity        cigar   comments
3ea026df-19e1-4102-9a6d-edb01697fc87    22.7    1.6     105     116     TCTTCGGGGG      TTTTTTTTTTT     78      1.0     39=3S
e6294bc0-b45b-459e-bf14-25c9eac0f6bf    24.1    1.9     110     120     TCTTCGGGGG      TTTTTTTTTT      73      0.974   21=1I17=3S
3bb32870-591a-419f-a2f8-fea85f7117fc    30.3    2.0     104     130     TCTTCGGGGG      TTTTTTTTTTTTTTTTTTTTTTTTTT      78      1.0     39=3S
13a240f9-71f8-4178-a4b0-d53a2dc8b36b    27.7    2.1     115     128     TCTTCGGGGG      TTTTCTTTTTTTT   55      0.872   22=1D3=1X1=1I1X1=1X8=3S
0da31119-fbef-4009-98b5-9303f6673b73    21.5    1.9     109     133     TCTTCGGGGG      TTTTTTTTTTTTTTTTTTTTTTTT        78      1.0     39=3S
```

- 30 bases long poly-A tails with 5 terminal G bases (-GGGGG)
```bash
read_id pt_length       per_base        pt_start        pt_end  before_pt       pt_seq  score   identity        cigar   comments
06780cfc-a5b2-4401-8835-84fa5a355378    22.9    2.1     110     118     ATCTTCCCCC      TTTTTTTT        78      1.0     39=3S
6828e8ac-b3c2-42a3-b6a7-205a3d62616c    31.1    1.9     107     115     TCTTCCCCCC      TTTTTTTT        56      0.872   11=2I12=3I11=3S
7eb56aa6-1467-4e0c-9539-767d76112d3d    23.8    2.1     111     128     CTCCCCCCCC      TTTTTTTTTTTTTTTTT       58      0.892   13=1X3=1D1=1I13=1D5=5S
cb70c6c7-e9f8-45a1-ad55-a3a3937a4a09    31.3    1.7     99      103     TCGCTCTACC      TTTT    72      0.951   34=1X3=1I2=1S
d8e413ad-1587-4be6-8b3e-8827a9ff3a76    33.1    1.8     100     109     TCTTCCCCCC      TTTTTTTTT       75      0.974   15=1D24=3S
```

2. New N3Pseq oligo `CAGCACCT ACTTGCCTGTCGCTCTATCTGCAGAGCAGAG TTT`

- yeast total RNA
```bash
read_id pt_length       per_base        pt_start        pt_end  before_pt       pt_seq  score   identity        cigar   comments
1a822816-2cff-4bfe-b622-be142ada2eb8    -10.0   2.2     -1      0       TTGCCGACTT              73      0.974   27=1I11=3S      25s     rRNA
ece6abde-9f6d-4680-95a8-82a64603d368    1.0     1.6     124     126     GTAATGATCC      TT      67      0.9     13=1D3=1X2=2D21=2S      18s     rRNA
29964a2e-80e4-4800-bfc5-0de8d724ae4d    20.1    1.9     98      115     CAGAGCAGAG      TTTTTTTTTTTTTTTTT       84      1.0     42=     YDR002W protein_coding
03d98db8-d6c9-459c-92c5-eeb0bec5f7fe    -10.0   2.1     -1      0       CTGCTTCGGT              78      1.0     39=3S   25s     rRNA
a89995a7-4844-48b9-a6e2-84b700745580    -2.4    2.4     106     108     AGCAGAGTAA      TT      75      0.975   11=1X28=2S      SCR1    ncRNA
```



## Citation

