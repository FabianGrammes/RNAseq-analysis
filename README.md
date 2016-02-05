# RNA_map

`RNA_map` is a bash script that automates the processes of:

1. Trimming of adapter sequences and poor quality nucleotides
   (cutadapt); Removing read pairs where one of the reads has less than 40bp 
2. Quality Control (FastQC)
3. Mapping reads to the genome (STAR)
4. Read counting per gene (HTSeq)
5. Summary statistics
6. Copy to common folder

Table of Contents
=================

  * [RNA_map](#rna_map)
    * [Example](#example)
      * [Reports](#reports)
    * [RNA_map options](#rna_map-options)
    * [<em>master sheet</em> format guide](#master-sheet-format-guide)
      * [TruSeq adaptors](#truseq-adaptors)
        * [Example 1: Same sample different lanes](#example-1-same-sample-different-lanes)
        * [Example 2: Single adaptors](#example-2-single-adaptors)
    * [<em>description file</em> format guide](#description-file-format-guide)
  * [Changelog](#changelog)
    * [Version 1.0.1](#version-101)
    * [Version 1.0](#version-10)

Created by [gh-md-toc](https://github.com/ekalinin/github-markdown-toc)

The script currently works only if the following conditions are met
(possible to extend the script for other settings):
- Libraries were prepared using a *TruSeq stranded Kit*
- fastq data was generated on a MiSeq/NextSeq/HiSeq

## Example

It's usually a *good idea* to run the script with the option: `--execute
no`; Check that the *number of samples* as well as the ID of the first
and last sample are correct; and than run the whole stuff. 

```bash
script=/mnt/users/fabig/cluster_pipelines/RnaMapping/RNA_map.sh
fdir=/mnt/SeqData2/MiSeq_data/150612_M02210_0008_000000000-ACGHJ/Data/Intensities/BaseCalls
sdir=/mnt/SeqData2/MiSeq_data/150612_M02210_0008_000000000-ACGHJ/Data/Intensities/BaseCalls/SampleSheet.csv

bash $script -d $fdir -s $sdir -m Samples.txt -r long --copy test --desc description.txt --execute no
```

### Reports

The script also submits 2 `R` scripts that generate both:

1. STAR mapping statistics (see: `star/STAR_Log_Stat.txt` and `star/STAR_Log_Stat.pdf`)
2. HTseq counting statistics (see: `count/HTSeq_Count_Stats.txt` and
`count/HTSeq_Count_Stats.pdf`)
3. Summarising read quality for untrimmed and trimmed reads

In order for the `R` scripts to work you need to have the following
packages installed in your `orion R`

- `ggplot2`
- `reshape`
- `gtable`
- `plyr`
- `grid`
- `gridExtra`

## RNA_map options

- `-d|--dirin`: *Required* Full path of the directory that holds the .fasq.gz
  files

- `-g|--genome`: *Required*, Path to STAR index of the Genome to be used

- `--gtf`: *Required*, Path to the .gtf/.gff file corresponding to the
  Genome. The .gtf/.gff information is used for counting the reads per
  genes. 

- `-m|--mastersheet`: *Required*  needs to point to a _master sheet_
file (see [ _master sheet_ format guide])

- `--idx`:*Optional* `<single/dual>` Defaults to `<single>`. Weather
single or dual indexes  have been used. Currently works only with
TruSeq aadptors see [ _master sheet_ format guide] on how to set up a
master sheet for single/dual indexed reads. 

- `-r|--read`: *Optional*  `<short/long>` Options passed on to STAR.
  Use <long> if the reads are longer than 2x250bp. Defaults to `<short>`.

- `--stranded`:*Optional*  `<yes/no/reverse>` Options passed on to
  `htseq-count` script. Defaults to `<reverse>` (for illumina stranded
  data); see
  [HtSeq](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html).

- `--trimmer`: *Optional*  `<cutadapt/trimmomatic>`; Defaults to `<cutadapt>`

- `--copy`: *Required* Name of the common folder where a copy of
  - count
  - mapping_summary
  - mastersheet
  will be placed. When set to `<no>`, nothing will be copied

- `--desc`: *Required* Short description file giving the basic informations about
  the experiment ( see [format guide](## _description file_ format guid)). **NOTE** If you
  use the `--copy no` then you do *NOT* need to provide file here. 

- `--execute`: If set to `<no>` all folders/scripts will be genearted,
  but none of the jobs will be executed. 

---

`RNA_map` creates a folder tree in the current directory:

- `slurm` folder to collect the slurm reports. *Note* each slurm report lists:
   - date of job submission
   - version of the used module 
- `bash` folder to collect the bash scripts, generated and submitted through `RNA_map`.
- `fastq_trim` folder for the adaptor and quality trimmed .fastq
  file, read pairs where  one 1 of the reads
  was found to be too short ( < 40bp ) where removed. Necessary for
  STAR to function properly. 
- `qc` folder for the quality control reports of the .fastq files.
- `qc_trim` folder for the quality control reports of the quality trimmed .fastq files.
- `star` folder contains read alignments to the genome in `.BAM`
format. Reads were aligned using the program `STAR`
- `count` folder contains read counts per gene generated by the `HTSeq-count` script. 
- `mapp_summary` folder that contains summarie statistics from STAR
  and HTSeq.

## _master sheet_ format guide

The script currently only works with Illumina TruSeq adaptors both single
or dual index.


### TruSeq adaptors

Plain tab delimited text file, no quotes, include header. *Important*
writing your _master sheet_ in excel works. If you have several .fastq
files per Sample, say comming from sequerncing the same sample on
different HiSeq lanes, The script will quality trim the .fastq files
individual and then summarise them into one .bam file (and
subsequently produce one .count file sample). In order for this to
work all files belonging to same sample need to have the same base
name: 

Using single or dual adaptors (TruSeq) works. See
[Illumina adaptor](https://www.med.unc.edu/pharm/calabreselab/files/tufts-sequencing-primer)
for further information about TruSeq adaptors.

#### Example 1: Same sample different lanes

- Sample 1; sequenced on Lane 1 and 3: Adaptor p7 = A001 (ATCACG)
	- `0-1_S1_L001_R1.fastq.gz`
	- `0-1_S1_L001_R2.fastq.gz`
   	- `0-1_S1_L003_R1.fastq.gz`
	- `0-1_S1_L003_R2.fastq.gz`
- Sample 2; sequenced on Lane 2 and 5: Adaptor A002 (CGATGT)
	- `CTR-1_S2_L002_R1.fastq.gz`
	- `CTR-1_S2_L002_R2.fastq.gz`
	- `CTR-1_S2_L005_R1.fastq.gz`
	- `CTR-1_S2_L005_R2.fastq.gz`

The _master sheet_ has to look like this (*WITH header*), Important it
is NOT necessary to include adaptor information in the _master_sheet_;
however it is good idea to do so anyway. **I recommend to use on column
per adaptor**. 

sample |	base 
----|-----
0_1 |  0-1_S1
CTR_1 |  CTR-1

#### Example 2: Single adaptors

Given 3 samples with single adaptors

- Sample 1: Adaptor p7 = A001 (ATCACG)
	- `0-1_S1_L001_R1_001.fastq.gz`
	- `0-1_S1_L001_R2_001.fastq.gz`
- Sample 2: Adaptor A002 (CGATGT)
	- `CTR-1_S2_L001_R1_001.fastq.gz`
	- `CTR-1_S2_L00_R2_001.fastq.gz`
- Sample 3: Adaptor A003 (TTAGGC)
	- `MA-1_S3_L001_R1_001.fastq.gz`
	- `MA-1_S3_L001_R2_001.fastq.gz`

The _master sheet_ has to look like this (*WITH header*), Important it
is NOT necessary to include adaptor information in the _master_sheet_;
however it is good idea to do so anyway. **I recommend to use on column
per adaptor**. 

sample |	base 
----|-----
0_1 |  0-1_S1_L001 
CTR_1 |  CTR-1_S2_L001
MA_1  | MA-1_S3_L001


Recomendet form:

sample |	base |	p7.adapt
----|-----|-----
0_1 |  0-1_S1_L001 |	A001 
CTR_1 |  CTR-1_S2_L001 | A002 
MA_1  | MA-1_S3_L001|  A003 



## _description file_ format guide

The following information should be contained in the description file:

1. Organism
2. Tissue
3. Date, contact person
4. 2-3 Sentences about goal/setup of the experiment

# Changelog

## Version 1.0.1

- Accepts multiple .fastq files per Sample as input. Several .fastq
  files will be summarised into one .bam file. All .fastq files need
  to have the same file prefix. 
  
- *STAR*
	- Not running as loop but now as array job with 5CPUs per task and
      max 20 jobs in parallel.
- *FASTQ*
	- Added R script to summarise the FastQC reports. Both .pdf and .txt
  
## Version 1.0

- Simplified the design of the master_sheet
- Modified the cutadapt script to work with all kinds of Illumina
  adators (single/dual adaptors). 
