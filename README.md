## Updating RepEnrich2 for Python 3.8 and STAR (RepEnrich2.38)

A complete rewrite of [RepEnrich2](https://github.com/nerettilab/RepEnrich2) which is a continuation of [RepEnrich](https://github.com/nskvir/RepEnrich) based on [this publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4122776/). **I do not take any credit for the methods**, this repository is an extension of their package. I found the original source code to be difficult to modify and thus difficult to directly implement new features in, so instead of a fork, I decided to rewrite the package from scratch. `fastRepEnrich` is written in Python 3.8 and utilizes the STAR aligner. The main speed increases come from using [STAR](https://github.com/alexdobin/STAR) instead of [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), utilizing Python builtins more effectively, and implementing multiprocessing during the pseudogenome building step.

## TODO:

- Benchmark versions
- Compare RepEnrich2 to fastRepEnrich results on real data
- add tests

## Installation

To ensure that all dependencies are met, first install and activate the conda environment using the environment.yml file.

`conda env create -f environment.yml`

This will create a conda environment called `fastre`. Activate the environment.

`conda activate fastre`

## Basic Workflow for Paired End Reads

Below I will detail my workflow with fastRepEnrich, using chr22 as an example and some simulated paired-end reads. The workflow for single end reads is basically the same, simply omit the second read and drop the `--pairedEnd` flags in fastRE commands.

Create a project directory with the following structure:

```
project/
	fastRE_count.py
	fastRE_setup.py
	fastRE_subset.py
	data/
		sample1/
			sample1_R1.fastq
			sample1_R2.fastq
		sample2/
			sample2_R1.fastq
			sample2_R2.fastq
		sampleN/
			...
		chr22.fa               # genome fasta file
		chr22.fa.out           # repeat masker file
		chr22.knownGenes.gtf   # annotation file
		sample-names.txt
```

Where sample-names.txt contains:

```
sample1
sample2
sample3
sample4
sample5
sample6
```

## 1. Create STAR idx for genome of interest

This will be tailored to your reads. **Since we need to prohibit splicing, you should build the genome index without using a GTF file**. See the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) for more details.

```bash
STAR --runThreadN 8 \
	--runMode genomeGenerate \
	--genomeDir data/STAR_chr22_idx \
	--genomeFastaFiles data/chr22.fa \
	--genomeSAindexNbases 11             # special setting for chr22
```

## 2. Map your samples to the genome, saving individual alignment results in each sample folder

Using STAR for this kind of mapping does not confer any additional benefits over bowtie2 aside from speed. Please see the following links for more information. [link1](https://groups.google.com/forum/#!searchin/rna-star/CHIP-seq%7Csort:date/rna-star/Gq3Gf3NmDNc/yt8uJ1u4AQAJ) [link2](https://groups.google.com/forum/#!searchin/rna-star/repeat$20elements%7Csort:date/rna-star/TqOdXiEFYrI/tbzoK_AV4DQJ) [link3](https://www.biostars.org/p/344389/)

The `--outFilterMultimapNmax 100` flag allows reads to align to N (100) different loci before being flagged as unmapped. **This flag is necessary for downstream analysis as it allows for multimapping reads (the default is set to 20 if not explicitly set)**. Additionally, in order to mimic the behavior of `bowtie2` we set the alignment type to end-to-end wth the `--alignEndsType EndToEnd` flag. The default `STAR` alignment method is local with soft-clipping allowed. We also set the `--alignIntronMin 1` flag to prohibit splicing.

```bash
for samp in $(cat data/sample-names.txt); do
	STAR --runThreadN 8 \
	--genomeDir data/STAR_chr22_idx/ \
	--readFilesIn data/${samp}/${samp}_R1.fastq data/${samp}/${samp}_R2.fastq \
	--outSAMtype BAM Unsorted \
	--outFilterMultimapNmax 100 \
	--outMultimapperOrder Random \
	--outFileNamePrefix data/${samp}/${samp}_ \
	--alignEndsType EndToEnd
	--alignIntronMin 1;
done
```

## 3. Run the fastRE_setup script

`python fastRE_setup.py data/chr22.fa.out data/chr22.fa --threads 8`

Here, you can use either the pre-assembled repeatmasker file downloaded from [repeatmasker.org](http://repeatmasker.org/) as is used above, or a custom bed file with the following format:

`chromosome\tstart_position\tend_position\trepeat_name\trepeat_class\trepeat_family`

```
chr22	10510228	10510528	AluSx1	SINE	Alu
chr22	10512707	10514778	L1MB1	LINE	L1
chr22	10514779	10515050	AluSx1	SINE	Alu
chr22	10515051	10515074	L1MB1	LINE	L1
chr22	10515075	10515121	(GAAG)n	Simple_repeat	Simple_repeat
chr22	10515122	10516103	L1MB1	LINE	L1
chr22	10516115	10516222	(TA)n	Simple_repeat	Simple_repeat
chr22	10516224	10516285	LTR66	LTR	ERVL
chr22	10516288	10516630	L1MB1	LINE	L1
chr22	10516636	10517247	L2a	LINE	L2
chr22	10517291	10517437	L1MEh	LINE	L1
chr22	10517525	10517865	L1M6	LINE	L1
chr22	10518784	10519114	MLT1A0	LTR	ERVL-MaLR
chr22	10519674	10519746	FRAM	SINE	Alu
chr22	10519747	10519816	MER52A	LTR	ERV1
```

If using the custom bed file then the `--isBed` flag should be set in you fastRE_setup.py command like so:

`python fastRE_setup.py data/chr22_repeats.bed data/chr22.fa --threads 8 --isBed`

## 4. Run the fastRE_subset script on all samples

**Importantly, the `--threads` flag in the subset command refers to the number of additional threads to use, not the total number of threads.** So if you have 8 threads total then set this number to 7. The default value is set to 0.

```bash
for samp in $(cat data/sample-names.txt); do
	python fastRE_subset.py data/${samp}/${samp}_Aligned.out.bam ${samp} --threads 7 --pairedEnd;
done
```

## 5. Run fastRE_count.py on all samples to get the count data

The `--summarize` flag below includes counts at both the class and family levels for each counting strategy in addition to counts at the repeat_name level for each counting strategy (the default output). 

```bash
for samp in $(cat data/sample-names.txt); do
	python fastRE_count.py ${samp} fastRE_Setup/ data/${samp}/${samp}_unique.bam --pairedEnd --threads 8 --summarize;
done
```

Each sample sub-folder will now contain all of the count data produced by fastRE as .tsv files. Use this data as input for differential expression analysis. The sample folders will also contain all of the information from the alignment steps contained in the `sampleName_Log.final.out` files.