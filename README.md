# Reprogrammed EMIRGE
Reprogrammed EMIRGE is a method for reconstructing the true sequences in an environment based on metagenomic sequencing reads, using an expectation-maximization algorithm. 
Building upon EMIRGE (https://github.com/csmiller/EMIRGE), Reprogrammed EMIRGE has several improvements. It now supports user-defined reference databases and has been demonstrated to effectively construct full-length sequences of microbial functional genes in Haima cold seep.
If you use Refined EMIRGE in your work, please cite these manuscripts as appropriate.
> Phylogenetic diversity of functional genes in the deep-sea cold seep: a novel perspective on metagenomics (to be published).
> Miller CS, Baker BJ, Thomas BC, Singer SW, Banfield JF (2011) EMIRGE: reconstruction of full-length ribosomal genes from microbial community short read sequencing data. Genome biology 12: R44. doi:10.1186/gb-2011-12-5-r44.
> Miller CS, Handley KM, Wrighton KC, Frischkorn KR, Thomas BC, Banfield JF (2013) Short-Read Assembly of Full-Length 16S Amplicons Reveals Bacterial Diversity in Subsurface Sediments. PloS one 8: e56018. doi:10.1371/journal.pone.0056018.

## How to Start
Users need to download and install EMIRGE, as well as the necessary tools and packages required for running EMIRGE.
#### ○ python3
#### ○ python (tested with version 2.6), with the following packages installed:
  -BioPython
  -Cython
  -pysam
  -scipy / numpy
 #### ○ usearch
 #### ○ samtools
 #### ○ bowtie2 (required by Reprogrammed EMIRGE)
 #### ○ bowtie (required if you use EMIRGE)
 #### ○ vsearch

## Make the Reference Database
The `make_my_db.py` script is used for constructing the initial reference database, replacing emirge_makedb.py in EMIRGE.
Run the following for help:
`python3 make_my_db.py --help`


## Run Reprogrammed EMIRGE
The `remirge.py` script is the main program for Reprogrammed EMIRGE. Prior to running it, please manually set the absolute path of **sam_filter.py** in four locations.
For example, change 
`python3 sam_filter.py /dev/stdin /dev/stdout` 
to 
`python3 /home/usr/tools/EMIRGE/sam_filter.py /dev/stdin /dev/stdout`.
Run the following for help:
`python3 make_my_db.py --help`

**Note that the interpreter for `remirge.py` is python2, while the interpreter for all other python scripts are python3.**

`abundance_calculate.py` **and** `normalization.py` **are used for subsequent statistical analysis of the output files from** `remirge.py`**, as demonstrated in the article** *Phylogenetic diversity of functional genes in the deep-sea cold seep: a novel perspective on metagenomics.* **The scripts are provided here for reference.**
