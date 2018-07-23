---
layout: "default"
title: Setup
weight : 2
#permalink: /Setup/
---

<h1><p style="text-align: center">Setup</p></h1>

-----
<br>

vSNP INSTALLATION
=================

## Python environment setup

Script is written in Python and must be ran using Python 3.  Currently tested with Python 3.6. 

Anaconda is a highly trusted Python package distrubution platform.  

If Anaconda is not yet installed follow the Anaconda instructions to install on your platform.

    https://www.anaconda.com/download/     
    
When installing from the command line use Anaconda's default installation except... when asked to, "prepend to PATH", choose yes.
    
Once Anaconda is installed close and reopen your terminal.

---

<strong>Optional:</strong>
If you are currently using Python 2, and wish to keep it as your default, a virtual environment can be built.  `vSNP.py` can run in this virtual environment.  See notes at the bottom of this page.  It is recommend to not run in a virtual environment if possible.

---

<br>
Setup Bioconda channels.  Add them in the order shown below.  Order is important.

    ~$ conda config --add channels conda-forge
    ~$ conda config --add channels defaults # Can ignored warning
    ~$ conda config --add channels r
    ~$ conda config --add channels bioconda
    
Install programs.
    
    ~$ conda install pyvcf biopython bwa samtools picard abyss gatk raxml newick_utils xlrd xlsxwriter gitpython regex pandas

Install cairosvg.  At the time of writing this (2018-02-22) cairosvg via conda was not compatiable with Python 3.6.

    ~$ pip install cairosvg

Cairosvg requires cairo.  If importing cairosvg fails, it is needed.  Follow instructions at: https://www.cairographics.org/download/

When gatk is downloaded using Anacoda it still needs to be registered.  GATK has a way to do this.  Go to GATK's website, download the GATK package, unzip it, and run:

    ~$ gatk-register /path/to/Downloads/GenomeAnalysisTK*/GenomeAnalysisTK.jar
    
After `gatk-register` is ran, GATK just downloaded from the GATK website, can be deleted.  The download was only needed to register the Anaconda GATK package.

## Script and dependents
Clone script to home directory: 

    ~$ git clone https://github.com/USDA-VS/snp_analysis.git
    
If git is unavailable, `~$ conda install git` will make it available.

Put `snp_analysis` in your $PATH, or easier run lines below to put script in your anaconda PATH.

    ~$ ln -s ~/snp_analysis/vSNP.py ~/anaconda*/bin/

## Dependency files
vSNP requires reference files, high quality VCF files for GATK best practices, defining SNP positions to group closely related isolates and those positions to filter.  These files are provide as dependency files.  When vSNP is ran it looks for these file in your home directory.  If they are not found vSNP will attempt to download the files from GitHub.  These files can also be added now.

    ~$ cd ~; git clone https://github.com/USDA-VS/dependencies.git

## Step 1 test

Test files can be downloaded at:

    ~$ git clone https://github.com/USDA-VS/fastq_data_set-tb_complex
    ~$ git clone https://github.com/USDA-VS/fastq_data_set-brucella
    
Files have been cut to 200,000 reads, which will give around 20X coverage.  This file size is convenent for downloading and testing.  They should not be added to any currated database or used in reporting.  The complete sequence files are available in SRA.

Test by making directory containing FASTQs the working directory.

    ~$ cd ~/fastq_data_set-tb_complex

To aid in testing make a backup of the files

    $ mkdir original test; cp *gz original; mv *gz test; cd test; ls

`vSNP.py` must only see `*fastq.gz` files in the working directory.  It will exit if any other file type is found.  `vSNP` will batch FASTQs based on available computer resources.

    $ vSNP.py

If an error occurs it may have to do with running multiple samples and downloading dependencies.  Reset `test` directory and restart.

    $ rm -r ../test/*; cp ../original/* ./; ls
    $ vSNP.py

## Step 2 test

In step 1 vSNP saw FASTQ files in the working directory and ran appropriately.  In step 2 it looks for VCF files.  `vSNP.py` must only see `*vcf` files in the working directory.  It will exit if any other file type is found.  

Test using VCF test files, or better yet use the VCFs you just produced from step 1 above.  From script 1 output use VCF files in the `alignment` directory ending in *_zc.vcf.  Make a working directory containing only those VCFs and call `vSNP.py`.  
    
If using VCF test files

    $ cd ~
    
    ~$ git clone https://github.com/USDA-VS/vcf_test_files.git
    
    ~$ cd vcf_test_files/bovis

    $ vSNP.py
    
For list of options:
    
    $ vSNP.py -h
    
<br>

---

<br>

<strong>Optional:</strong>
If you are currently using Python 2, and wish to keep it as your default, a virtual environment can be built.  `vSNP.py` can run in this virtual environment.

To setup virtual environment:

    ~$ conda create -n snp_analysis python=3.6 # Optional:  skip if wanting to use your current environment
    ~$ source activate snp_analysis # Optional:  skip if wanting to use your current environment

Proceed to setting up channels...

To deactivate virtual environment, use:
    
    > source deactivate <env name> # skip if new environment was not activated
    
