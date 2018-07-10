#!/usr/bin/env python

import zipfile
import xlsxwriter
import xlrd
import vcf
import time
import sys
import subprocess
import smtplib,ssl
import shutil
import regex
import re
import numpy as np
import pandas as pd
import os
import multiprocessing
from multiprocessing import Pool
import gzip
import glob
import git
import csv
import argparse
import textwrap
import signal
import json
from collections import defaultdict
from collections import Iterable
from cairosvg import svg2pdf
from numpy import mean
from functools import partial
from email.utils import formatdate
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders
from distutils.dir_util import copy_tree
from datetime import datetime
from concurrent import futures
from collections import OrderedDict
from collections import Counter
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO


###############################################
###############################################
##################loopwrapper##################
###############################################
###############################################

class loop():
    
    def run_loop(self):
        home = os.path.expanduser("~")

        startTime = datetime.now()
        print ("\n\n*** START ***\n")
        print ("Start time: %s" % startTime)

        list_of_files = glob.glob('*gz')
        list_len = len(list_of_files)
        if (list_len % 2 != 0):
            print("\n#####Check paired files.  Unpaired files seen by odd number of counted FASTQs\n\n")
            sys.exit(0)

        for file in list_of_files:
            prefix_name=re.sub('_.*', '', file)
            print(prefix_name)
            if not os.path.exists(prefix_name):
                os.makedirs(prefix_name)
            shutil.move(file, prefix_name)

        ###
        #Run stats

        ts = time.time()
        st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')

        #placed at root
        #get file opened and give a header
        summary_file = root_dir + '/stat_alignment_summary_' + st + '.xlsx'
        workbook = xlsxwriter.Workbook(summary_file)
        worksheet = workbook.add_worksheet()
        row = 0
        col = 0
        top_row_header = ["time_stamp", "sample_name", "self.species", "reference_sequence_name", "R1size", "R2size", "Q_ave_R1", "Q_ave_R2", "Q30_R1", "Q30_R2",  "allbam_mapped_reads", "genome_coverage", "ave_coverage", "ave_read_length", "unmapped_reads", "unmapped_assembled_contigs", "good_snp_count", "mlst_type", "octalcode", "sbcode", "hexadecimal_code", "binarycode"]
        for header in top_row_header:
            worksheet.write(row, col, header)
            col += 1
        ###

        #Cumulative stats
        path_found = False
        if os.path.isdir("/bioinfo11/TStuber/Results/stats"): #check bioinfo from server
            path_found = True
            copy_to = "/bioinfo11/TStuber/Results/stats"
        elif os.path.isdir("/Volumes/root/TStuber/Results"): #check bioinfo from Mac
            path_found = True
            copy_to = "/Volumes/root/TStuber/Results/stats"
        else:
            copy_to="no_path"
            print("Bioinfo not connected")

        if path_found:
            try:
                summary_cumulative_file = copy_to + '/stat_alignment_culmulative_summary' + '.xlsx'
                summary_cumulative_file_temp = copy_to + '/stat_alignment_culmulative_summary-' + st + '-temp.xlsx'
                temp_folder = copy_to + '/temp'
            except OSError:
                print("\n\nBioinfo unresponsive\n\nUnable to copy to stats file\n\n")
                text = "ERROR, Bioinfo unresponsive unable to copy to stats file"
                msg = MIMEMultipart()
                msg['From'] = "tod.p.stuber@aphis.usda.gov"
                msg['To'] = "tod.p.stuber@aphis.usda.gov"
                msg['Date'] = formatdate(localtime = True)
                msg['Subject'] = "### No coverage file"
                msg.attach(MIMEText(text))
                smtp = smtplib.SMTP('10.10.8.12')
                smtp.send_message(msg)
                smtp.quit()
        ###

        directory_list=[]
        for f in  os.listdir('.'):
            if not f.startswith('.'):
                directory_list.append(f)

        total_samples = len(directory_list)
        lower_count = 0
        upper_count = 1
        row = 1
        while lower_count < total_samples:
            upper_count = lower_count + limited_cpu_count
            run_list = directory_list[lower_count:upper_count] #create a run list
            for i in run_list:
                directory_list.remove(i)
            total_samples = len(directory_list)
            print(run_list)

            print("Iterating directories")
            frames = []
            if args.debug_call: #run just one sample at a time to debug
                for d in run_list:
                    print("DEBUGGING, SAMPLES RAN INDIVIDUALLY")
                    stat_summary = read_aligner(d)
                    df_stat_summary = pd.DataFrame.from_dict(stat_summary, orient='index') #convert stat_summary to df
                    frames.append(df_stat_summary) #frames to concatenate

                    worksheet.write(row, 0, stat_summary.get('time_stamp', 'n/a'))
                    worksheet.write(row, 1, stat_summary.get('sample_name', 'n/a'))
                    worksheet.write(row, 2, stat_summary.get('self.species', 'n/a'))
                    worksheet.write(row, 3, stat_summary.get('reference_sequence_name', 'n/a'))
                    worksheet.write(row, 4, stat_summary.get('R1size', 'n/a'))
                    worksheet.write(row, 5, stat_summary.get('R2size', 'n/a'))
                    worksheet.write(row, 6, stat_summary.get('Q_ave_R1', 'n/a'))
                    worksheet.write(row, 7, stat_summary.get('Q_ave_R2', 'n/a'))
                    worksheet.write(row, 8, stat_summary.get('Q30_R1', 'n/a'))
                    worksheet.write(row, 9, stat_summary.get('Q30_R2', 'n/a'))
                    worksheet.write(row, 10, stat_summary.get('allbam_mapped_reads', 'n/a'))
                    worksheet.write(row, 11, stat_summary.get('genome_coverage', 'n/a'))
                    worksheet.write(row, 12, stat_summary.get('ave_coverage', 'n/a'))
                    worksheet.write(row, 13, stat_summary.get('ave_read_length', 'n/a'))
                    worksheet.write(row, 14, stat_summary.get('unmapped_reads', 'n/a'))
                    worksheet.write(row, 15, stat_summary.get('unmapped_assembled_contigs', 'n/a'))
                    worksheet.write(row, 16, stat_summary.get('good_snp_count', 'n/a'))
                    worksheet.write(row, 17, stat_summary.get('mlst_type', 'n/a'))
                    worksheet.write(row, 18, stat_summary.get('octalcode', 'n/a'))
                    worksheet.write(row, 19, stat_summary.get('sbcode', 'n/a'))
                    worksheet.write(row, 20, stat_summary.get('hexadecimal_code', 'n/a'))
                    worksheet.write(row, 21, stat_summary.get('binarycode', 'n/a'))
                    row += 1

                    os.chdir(root_dir)
            else: # run all in run_list in parallel
                print("SAMPLES RAN IN PARALLEL")
                with futures.ProcessPoolExecutor(max_workers=limited_cpu_count) as pool: #max_workers=cpu_count
                    for stat_summary in pool.map(read_aligner, run_list): #run in parallel run_list in read_aligner (script1)
                        df_stat_summary = pd.DataFrame.from_dict(stat_summary, orient='index') #convert stat_summary to df
                        frames.append(df_stat_summary) #frames to concatenate

                        worksheet.write(row, 0, stat_summary.get('time_stamp', 'n/a'))
                        worksheet.write(row, 1, stat_summary.get('sample_name', 'n/a'))
                        worksheet.write(row, 2, stat_summary.get('self.species', 'n/a'))
                        worksheet.write(row, 3, stat_summary.get('reference_sequence_name', 'n/a'))
                        worksheet.write(row, 4, stat_summary.get('R1size', 'n/a'))
                        worksheet.write(row, 5, stat_summary.get('R2size', 'n/a'))
                        worksheet.write(row, 6, stat_summary.get('Q_ave_R1', 'n/a'))
                        worksheet.write(row, 7, stat_summary.get('Q_ave_R2', 'n/a'))
                        worksheet.write(row, 8, stat_summary.get('Q30_R1', 'n/a'))
                        worksheet.write(row, 9, stat_summary.get('Q30_R2', 'n/a'))
                        worksheet.write(row, 10, stat_summary.get('allbam_mapped_reads', 'n/a'))
                        worksheet.write(row, 11, stat_summary.get('genome_coverage', 'n/a'))
                        worksheet.write(row, 12, stat_summary.get('ave_coverage', 'n/a'))
                        worksheet.write(row, 13, stat_summary.get('ave_read_length', 'n/a'))
                        worksheet.write(row, 14, stat_summary.get('unmapped_reads', 'n/a'))
                        worksheet.write(row, 15, stat_summary.get('unmapped_assembled_contigs', 'n/a'))
                        worksheet.write(row, 16, stat_summary.get('good_snp_count', 'n/a'))
                        worksheet.write(row, 17, stat_summary.get('mlst_type', 'n/a'))
                        worksheet.write(row, 18, stat_summary.get('octalcode', 'n/a'))
                        worksheet.write(row, 19, stat_summary.get('sbcode', 'n/a'))
                        worksheet.write(row, 20, stat_summary.get('hexadecimal_code', 'n/a'))
                        worksheet.write(row, 21, stat_summary.get('binarycode', 'n/a'))
                        row += 1

                if not args.quiet and path_found:
                    try:
                        open_check = open(summary_cumulative_file, 'a') #'a' is very important, 'w' will leave you with an empty file
                        open_check.close()
                        df_all=pd.read_excel(summary_cumulative_file)
                        df_all_trans = df_all.T #indexed on column headers
                        # save back the old and remake the working stats file
                        shutil.move(summary_cumulative_file, '{}' .format(temp_folder + '/stat_backup' + st + '.xlsx'))
                        sorter = list(df_all_trans.index) #list of original column order
                        frames.insert(0, df_all_trans) #put as first item in list
                        df_concat = pd.concat(frames, axis=1) #cat frames
                        df_sorted = df_concat.loc[sorter] #sort based on sorter order
                        df_sorted.T.to_excel(summary_cumulative_file, index=False) #transpose before writing to excel, numerical index not needed
                    except BlockingIOError:
                        sorter = list(df_stat_summary.index) #list of original column order
                        df_concat = pd.concat(frames, axis=1) #cat frames
                        df_sorted = df_concat.loc[sorter] #sort based on sorter order
                        df_sorted.T.to_excel(summary_cumulative_file_temp, index=False)
                    except OSError:
                        sorter = list(df_stat_summary.index) #list of original column order
                        df_concat = pd.concat(frames, axis=1) #cat frames
                        df_sorted = df_concat.loc[sorter] #sort based on sorter order
                        try:
                            df_sorted.T.to_excel(summary_cumulative_file_temp, index=False)
                        except OSError:
                            print("##### UNABLE TO MAKE CONNECTION TO BIOINFO")
                            pass
                else:
                    print("Path to cumulative stat summary file not found")

####send email:
        def send_email(email_list, runtime):
            text = "See attached:  "
            send_from = "tod.p.stuber@aphis.usda.gov"
            send_to = email_list
            msg = MIMEMultipart()
            msg['From'] = send_from
            msg['To'] = send_to
            msg['Date'] = formatdate(localtime = True)
            if not path_found:
                msg['Subject'] = "###CUMULATIVE STATS NOT UPDATED - Script1 stats summary"
            else:
                msg['Subject'] = "Script1 stats summary, runtime: {}" .format(runtime)
            msg.attach(MIMEText(text))

            part = MIMEBase('application', "octet-stream")
            part.set_payload(open(summary_file, "rb").read())
            encoders.encode_base64(part)
            part.add_header('Content-Disposition', 'attachment; filename="summary_file.xlsx"')
            msg.attach(part)

            #context = ssl.SSLContext(ssl.PROTOCOL_SSLv3)
            #SSL connection only working on Python 3+
            smtp = smtplib.SMTP('10.10.8.12')

            smtp.send_message(msg)
            #smtp.sendmail(send_from, send_to, msg.as_string())
            smtp.quit()

        workbook.close()

        runtime = (datetime.now() - startTime)
        print ("\n\nruntime: %s:  \n" % runtime)

        if args.email:
            send_email(email_list, runtime)



################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
def get_species():

    #species = corresponding NCBI accession
    species_cross_reference = {}
    species_cross_reference["salmonella"] = ["016856, 016855"]
    species_cross_reference["bovis"] = ["AF2122_NC002945", "00879"]
    species_cross_reference["af"] = ["NC_002945.4"]
    species_cross_reference["h37"] = ["000962", "002755", "009525", "018143"]
    species_cross_reference["para"] = ["NC_002944"]
    species_cross_reference["ab1"] = ["006932", "006933"]
    species_cross_reference["ab3"] = ["007682", "007683"]
    species_cross_reference["canis"] = ["010103", "010104"]
    species_cross_reference["ceti1"] = ["Bceti1Cudo"]
    species_cross_reference["ceti2"] = ["022905", "022906"]
    species_cross_reference["mel1"] = ["003317", "003318"]
    species_cross_reference["mel1b"] = ["CP018508", "CP018509"]
    species_cross_reference["mel2"] = ["012441", "012442"]
    species_cross_reference["mel3"] = ["007760", "007761"]
    species_cross_reference["ovis"] = ["009504", "009505"]
    species_cross_reference["neo"] = ["KN046827"]
    species_cross_reference["suis1"] = ["017250", "017251"]
    species_cross_reference["suis2"] = ["NC_010169", "NC_010167"]
    species_cross_reference["suis3"] = ["007719", "007718"]
    species_cross_reference["suis4"] = ["B-REF-BS4-40"]
    
    vcf_list = glob.glob('*vcf')
    for each_vcf in vcf_list:
        print(each_vcf)
        mal = fix_vcf(each_vcf)
        vcf_reader = vcf.Reader(open(each_vcf, 'r'))
        print("single_vcf %s" % each_vcf)
        for record in vcf_reader:
            header = record.CHROM
            for k, vlist in species_cross_reference.items():
                for l in vlist:
                    if l in header:
                        return(k)

global root_dir
root_dir = str(os.getcwd())

global cpu_count
global limited_cpu_count
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
print("")

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

if fastq_check > 0:
    fastq_check = True
if vcf_check > 0:
    vcf_check = True
if fastq_check and vcf_check:
    print("\n#####You have a mix of FASTQ and VCF files.  This is not allowed\n\n")
    sys.exit(0)

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
            print("\n--> RUNNING LOOP/SCRIPT 1\n") #
            loop().run_loop()
elif vcf_check:

    if not args.species:
        args.species = get_species()
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
                print("\n--> RUNNING SCRIPT 2\n") #
                script2().run_script2()
            else:
                print("#####Based on VCF CHROM id (reference used to build VCF) a matching species cannot be found neither was there a -s option given")
                sys.exit(0)

else:
    print ("#####Error determining file type.")
    sys.exit(0)

