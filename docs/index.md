---
# You don't need to edit this file, it's empty on purpose.
# Edit theme's home layout instead if you wanna make some changes
# See: https://jekyllrb.com/docs/themes/#overriding-theme-defaults

title: Overview
layout: "default"
weight: 1
---

<h1><p style="text-align: center">Overview</p></h1>

-----
<br>
vSNP --> USDA APHIS Veterinary Services pipeline for Mycobacterium tuberculosis complex and Brucella sp.  Genotyping from high throughput sequence providing SNP tables and phylogentic trees with output to aid in SNP validation. 

vSNP can be installed on Linux and Mac platforms and is ran from the command-line.  It is written in Python 3 and relies on the Anaconda package manager.  It is used in a two step process.  Step 1 is called on a working directory containing FASTQ files.  BWA is used to align reads and SNPs are called using GATK's HaplotypeCaller outputing VCF files.  In step 2, those VCF files are gathered to produce SNP alignments, tables and phylogenetic trees of grouped isolates.  vSNP is portable and can be ran with relatively little computer resources.  Because vSNP groups samples of similar isolates one is able to quickly validate SNP positions, and report sample comparisons as a high-quality validated SNP analysis.

Minimal computer requirements are 4 cores, and 8GB of memory, but more compute resources are advantageous when running multiple samples, FASTQ file sizes are excessively large or there are over 1,000 VCFs in a comparison.

# Step 1 - FASTQ to VCF
Step 1 is fairly straight forward.  Our main workflow includes Mycobacterium tuberculosis complex and Brucella sp. and therefore the script has been optimized for such.  The script begins by selecting the best reference and determining the spoliogtype for TB complex isolates and MLST for Brucella spieces.  The script then goes on to perform what is now a fairly standard bacterial SNP calling pipeline where BWA aligns and GATK calls SNPs.  However, in addition to what is output by the HaplotypeCaller, zero coverage positions are added to the VCF files.  This is not a standard GATK option and is done using Python.  Including zero coverage positions allows a more accurate SNP summary to be represented in step 2.

# Step 2 - VCF to SNP alignment
Step 2 is called on VCF files output from step 1.  References chosen in step 1 have been selected because they have been found to be relatively close to the isolate.  The closer the reference is to the isolate the less overall SNP calling error is seen.  VCF files analyzed in step 2 must all be output from the same reference.  Obviously VCF files analyzed using different references can not be used in the same comparison.

In addition to using a closely related FASTA reference file, which minimizes SNP calling error, there are three additional external files, or dependencies, used to create high quality, informative SNP alignments.  As shown in bovis_dependency_view.jpg and suis1_dependency_view.jpg the reference used to build VCF files must be reflected in the three dependencies.  The three dependent files are: filter file, defining SNPs, and gbk file.  

To keep dependency files up-to-date in step 1 and 2, and to minimize the download requirement, files are downloaded to the user's home directory as needed when the script is ran.  The setup instructions also provides instructions for downloading files manually.