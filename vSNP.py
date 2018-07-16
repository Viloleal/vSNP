#!/usr/bin/env python

import os
import sys
import multiprocessing
import argparse
import textwrap
import glob

import functions

root_dir = str(os.getcwd())

#set cpu usage
cpu_count = multiprocessing.cpu_count()
limited_cpu_count = int(cpu_count/4)
if limited_cpu_count == 0:
    limited_cpu_count = 1

parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

---------------------------------------------------------
vSNP --> get SNPs, group SNPs, verify SNPs

vSNP is called on a working directory containing FASTQ or VCF files.

See documentation at: https://usda-vs.github.io/snp_analysis/

        Step 1: FASTQs --> VCF

        Step 2: VCFs --> Tables & Trees

-s <OPTIONAL SPECIES TYPES>: af, h37, ab1, ab3, suis1, suis2, suis3, mel1, mel1b, mel2, mel3, canis, ceti1, ceti2, ovis, neo, para, salmonella

'''), epilog='''---------------------------------------------------------''')

#universal
parser.add_argument('-s', '--species', action='store', dest='species', help='OPTIONAL: Used to FORCE species type <see options above>')
parser.add_argument('-d', '--debug', action='store_true', dest='debug_call', help='debug, run without loop.map for loops')
parser.add_argument('-g', '--get', action='store_true', dest='get', help='get, get to the core functions for debugging')
parser.add_argument('-n', '--no_annotation', action='store_true', dest='no_annotation', help='no_annotation, run without annotation')
parser.add_argument('-a', '--all_vcf', action='store_true', dest='all_vcf', help='make tree using all VCFs')
parser.add_argument('-e', '--elite', action='store_true', dest='elite', help='create a tree with on elite sample representation')
parser.add_argument('-f', '--filter', action='store_true', dest='filter', help='Find possible positions to filter')
parser.add_argument('-q', '--quiet', action='store_true', dest='quiet', help='[**APHIS only**] prevent stats going to cumlative collection')
parser.add_argument('-m', '--email', action='store', dest='email', help='[**APHIS only**, specify own SMTP address for functionality] email options: all, s, tod, jess, suelee, chris, email_address')
parser.add_argument('-u', '--upload', action='store_true', dest='upload', help='[**APHIS only**, specify own storage for functionality] upload files to the bioinfo drive')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.1')
args = parser.parse_args()
print ("\nSET ARGUMENTS: ")
print (args)
args_options = args
print("")

email_list = None
if args.email == "all":
    email_list = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Christine.R.Quance@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov, patrick.m.camp@aphis.usda.gov, David.T.Farrell@aphis.usda.gov, Robin.L.Swanson@aphis.usda.gov, hannah.m.tharp@aphis.usda.gov, Doris.M.Bravo@aphis.usda.gov, eto3@cdc.gov, kristina.lantz@aphis.usda.gov"
elif args.email == "tod":
    email_list = "tod.p.stuber@aphis.usda.gov"
elif args.email == "jess":
    email_list = "Jessica.A.Hicks@aphis.usda.gov"
elif args.email == "suelee":
    email_list = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov, Doris.M.Bravo@aphis.usda.gov, kristina.lantz@aphis.usda.gov"
elif args.email == "suelee-":
    email_list = "tod.p.stuber@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov"
elif args.email == "chris":
    email_list = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Christine.R.Quance@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov, eto3@cdc.gov, kristina.lantz@aphis.usda.gov"
elif args.email == "chris-":
    email_list = "tod.p.stuber@aphis.usda.gov, Christine.R.Quance@aphis.usda.gov"
elif args.email == "doris":
    email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, doris.m.bravo@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov, kristina.lantz@aphis.usda.gov"
else:
    email_list = "tod.p.stuber@aphis.usda.gov"

################################################################################################################################################

all_file_types_count = len(glob.glob('*.*'))
fastq_check = len(glob.glob('*fastq.gz'))
vcf_check = len(glob.glob('*vcf'))

# Check that there are either FASTQs or VCFs, not both
if fastq_check > 0:
    fastq_check = True
if vcf_check > 0:
    vcf_check = True
if fastq_check and vcf_check:
    print("\n#####You have a mix of FASTQ and VCF files.  This is not allowed\n\n")
    sys.exit(0)

# Check that there an equal number of both R1 and R2 reads
if fastq_check:
    R1 = glob.glob('*_R1*fastq.gz')
    R2 = glob.glob('*_R2*fastq.gz')
    R1count = len(R1)
    R2count = len(R2)
    fastq_count = R1count + R2count
    if (fastq_count % 2 != 0):
        print("\n#####Check paired files.  Unpaired files seen by odd number of counted FASTQs\n\n")
        sys.exit(0)
    if (R1count != R2count):
        print("\n#####Check paired files.  R1 files do not equal R2\n\n")
        sys.exit(0)
    if (all_file_types_count != fastq_count):
        print("\n#####Only zipped FASTQ files are allowed in directory\n\n")
        sys.exit(0)
    elif (fastq_count > 1):
        if args.all_vcf or args.elite or args.upload or args.filter:
            print("#####Incorrect use of options when running loop/script 1")
            sys.exit(0)
        else:
            print("\n--> RUNNING LOOP/SCRIPT 1\n")
            #Enter script 1 -->
            functions.run_loop(root_dir, limited_cpu_count, args_options, email_list)
            print("See files, vSNP has finished alignments")
elif vcf_check:
    if not args.species:
        args.species = functions.get_species()
        print("args.species %s" % args.species)
    vcfs_count = len(glob.glob('*vcf'))
    if (all_file_types_count != vcfs_count):
        print("\n#####You have more than just VCF files in your directory.  Only VCF files are allowed if running script 2\n\n")
        sys.exit(0)
    else:
        if args.quiet:
            print("#####Incorrect use of options when running script 2")
            sys.exit(0)
        else:
            if args.species:
                print("\n--> RUNNING SCRIPT 2\n")
                #Enter script 2 -->
                functions.run_script2()
            else:
                print("#####Based on VCF CHROM id (reference used to build VCF) a matching species cannot be found neither was there a -s option given")
                sys.exit(0)
else:
    print ("#####Error determining file type.")
    sys.exit(0)