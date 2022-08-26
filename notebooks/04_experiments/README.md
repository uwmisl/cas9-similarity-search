# This document explains how to analyze nanopore data for the Cas9 Similarity Search Project. 
# Setting up Docker is not required!

## General Workflow (AKA Quickstart)

First, you will need to establish what sequences you're examining by uploading a `.csv` file of all the sequences you're looking for. (More on how to do this below.)

Next, you will need to establish what nanopore data is your baseline data. You do not need to run this if you have already done so once. (More on how to do this below.) An example of this command is:

**YOU MUST GIVE THE FULL PATH FOR EACH ARGUMENT FOR THE SCRIPT TO WORK**


`python 1_css_initial_pool_analysis.py /Users/leeorg/Documents/MISL/cas9/nanopore_sim_search/seqs.csv /Users/leeorg/Documents/MISL/cas9/nanopore_sim_search/20220627_1431_MN40387_FAT53551_1a7f07f8/`

Finally, you will establish which nanopore run you wish to analyze and run a command in the command line that will take all this into account and output a PDF file with the results of various potentially useful analyses, as well as a directory of all the figures generated for easy use elsewhere. An example of this command is:

`python 2_css_run_analysis.py /Users/leeorg/Documents/MISL/cas9/nanopore_sim_search/seqs.csv /Users/leeorg/Documents/MISL/cas9/nanopore_sim_search/20220627_1431_MN40387_FAT53551_1a7f07f8/unmodified_pool.csv /Users/leeorg/Documents/MISL/cas9/nanopore_sim_search/11ultramers_20220808/20220808_1602_MN21390_FAT57906_9180b0f5 css_run_analysis.py`

## STEP 1: Setting up what sequences you want to use for alignment

### Why? 
So you know what reads to consider aligned to your targets. This file is used by both `css_pool_0_analysis.py` and `css_run_analysis.py`.

### How? 
Upload a `.csv` file where the first column is `seq_id` (just a name for each sequence, keep it short as these will be the names that appears on all figures) and the second column is `seq` (the DNA sequence). This file will be passed to the scripts it is used in, no need for a specific name, it just needs to be a `.csv` file with those two columns.

## STEP 2 (it's only necessary to do this step one time- the time you designate this run as a baseline run): Setting up which run is the initial run

### Why? 
So that your data can be normalized to the initial concentrations of each sequence.

### How?
Upload your nanopore data. The directory should include `fastq_pass` and `fastq_fail` directories (though it is not necessary to have both). These directories will include zipped `.fastq` files that will be unzipped and analyzed automatically, there is no need to unzip/modify them. 

To do this, you will run a command in the command line that takes the path to the nanopore data directory (one level up from the `fastq_pass` and `fastq_fail` directories) as well as the path to the `.csv` with sequences you're looking to align to.

An example of this command is in the quickstart section above.

### What is `1_css_initial_pool_analysis.py` doing?
This script takes in nanopore data and a pre-existing `.csv` file with sequences to align (to change this, you will need to change the global variable inside this python script). It outputs `unmodified_pool.csv` which contains every read that aligned to a sequence to align (read ID, highest alignment score, the read sequence, and which sequence it aligned to). 

## STEP 3: Analyzing subsequent run data

## How?
Upload your nanopore data. The directory should include `fastq_pass` and `fastq_fail` directories (though it is not necessary to have both). These directories will include zipped `.fastq` files that will be unzipped and analyzed automatically. 

Then run a command line command that takes the `unmodified_pool.csv` path and the path to your data as well as the path to the `.csv` with sequences you're looking to align to. An example of this command above in the quickstart section.


## What is `css_run_analysis.py` doing? 
First, all reads that align to a sequence of interest are identified. Then, all reads aligned to multiple sequences are discarded (the number of sequences discarded is noted in the final PDF). The number of sequences aligned are noted in the final PDF.

Next, the enrichment score for each sequence of interest is calculated (for more on the enrichment score, see Section 1 of the Supplemental in the published work). 

Finally, the sequence length distributions are examined.

This is all output in a final PDF with the format `yyyymmdd_css_analysis_summary.pdf` where the date comes from the nanopore directory name. All images are stored in a directory titled `yyyymmdd_css_analysis`. 
