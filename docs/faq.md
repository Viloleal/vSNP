---
layout: "default"
title: FAQ
weight: 7
permalink: /faq/
---

<h1><p style="text-align: center">Frequently Asked Questions</p></h1>

-----
<br>

## GATK fails

Did `gatk-register` run successfully?  See setup instructions.  After running command a message similar to "jar file specified matches expected version" should be output.

## Samtools unable to access libraries

`samtools: error while loading shared libraries: libbz2.so.1.0: cannot open shared object file: No such file or directory`

The Samtools Anaconda install is unable to access libraries.  Fix: Install from different conda-forge channel

`~$ conda install -c conda-forge -c bioconda samtools bzip2`


The above should now print samtools list of commands.

## Picard fails

Error:  Picard fails to find Java, causing script to fail.

The Anaconda version of Picard can be tested using:

`~$ picard`

This should output the Picard command line usage.  If an error occurs check Java version.  

`~$ java -version`

Both Picard and GATK require Java 1.8 (aka version 8).

Anaconda should have taken care of the your Java version.  There are different reasons which can cause Java to fail.  In my experience it has been caused by a previously set .bashrc (Linux) /.bash_profile (macOS) profile setting.  I recommend going into your profile and commenting out PATH settings that may no longer be needed.  Remember to restart your terminal to update profile befor testing again with: `~$ picard`

## Do I have to use Anaconda package manager?

No, the Anaconda package manager is not need, but it does make the setup much easier.  

System programs such as: bwa, samtools, abyss, and raxml will be usable if they are in your PATH with the exact same name (remember Linux is case sensitive meaning capitalization counts).  This is why samtools above can be removed from your conda environment and replaced with the standard install.  Either way `samtools` is being found in your PATH.  However Picard and GATK are Java programs.  If running from a standard installation the call within vSNP.py needs to be changed.  For example, when ran from a standard installation the call is: `java -Xmx1g -jar ${picard} CreateSequenceDictionary ...`.  The `java -Xmx1g -jar` must be included.  However, Anaconda takes care of this, which allows vSNP to use only `picard CreateSequenceDictionary ...` when Picard and GATK are install using it.

Also, the Anaconda package manager is not needed to install specific Python modules.  These can be installed using other tools such as pip or easy_install.  What is important is getting the correct version, i.e. Python 3.5 and pandas=0.18.1.

## How do I start over with the Anaconda installation?

Simply remove your anaconda folder, `rm -rf ~/anaconda`, close and reopen your terminal, and restart with the setup instructions. 

## Panda's error occurs at around line 3400.

The error may represent itself as, `KeyError: 'reference_pos'` or `*** ValueError: all the input arrays must have same number of dimensions`

This may be caused by either an old version of pandas or the latest version of the sript is not being used.

History:  Pandas 0.18.1 contained a bug which required a work around.  This bug has been fixed and the script has been updated properly append tables.

Fix:  Update to the latest verion of vSNP and pandas.  To update vSNP preform a git pull or re-clone the repo.  To update pandas: `~$ conda update pandas`.  Version must be > 0.18.1.  Tested with 0.20.3.

-----
