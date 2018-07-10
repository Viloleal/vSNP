#!/usr/bin/env python

import glob
import vcf
import re
import os
from concurrent import futures

def fix_vcf(each_vcf):
    mal = []
    ###
    # Fix common VCF errors
    # if args.debug_call:
        # print ("FIXING FILE: " + each_vcf)
    temp_file = each_vcf + ".temp"
    write_out=open(temp_file, 'w') #r+ used for reading and writing to the same file
    ###
    with open(each_vcf, 'r') as file:
        try:
            for line in file:
                if line.rstrip(): # true if not empty line'^$'
                    line = line.rstrip() #remove right white space
                    line = re.sub('"AC=', 'AC=', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('"$', '', line)
                    line = re.sub('GQ:PL\t"', 'GQ:PL\t', line)
                    line = re.sub('[0-9]+\tGT\t.\/.$', '999\tGT:AD:DP:GQ:PL\t1/1:0,80:80:99:2352,239,0', line)
                    line = re.sub('^"', '', line)
                    if line.startswith('##') and line.endswith('"'):
                        line = re.sub('"$', '', line)
                    if line.startswith('##'):
                        line = line.split('\t')
                        line = ''.join(line[0])
                    if not line.startswith('##'):
                        line = re.sub('"', '', line)
                        line = line.split('\t')
                        line = "\t".join(line[0:10])
                        print(line, file=write_out)
                    else:
                        print(line, file=write_out)
        except IndexError:
            print ("##### IndexError: Deleting corrupt VCF file: " + each_vcf)
            mal.append("##### IndexError: Deleting corrupt VCF file: " + each_vcf)
            os.remove(each_vcf)
        except UnicodeDecodeError:
            print ("##### UnicodeDecodeError: Deleting corrupt VCF file: " + each_vcf)
            mal.append("##### UnicodeDecodeError: Deleting corrupt VCF file: " + each_vcf)
            os.remove(each_vcf)

    write_out.close()
    os.rename(temp_file, each_vcf)
    return mal

vcf_list = glob.glob('*vcf')
print("Fixing files...\n")
#debug
malformed = []
for each_vcf in vcf_list:
    print(each_vcf)
    mal = fix_vcf(each_vcf)
    malformed = malformed + list(mal)
#pool
# with futures.ProcessPoolExecutor() as pool:
#     mal = pool.map(fix_vcf, vcf_list)
#     malformed = malformed + list(mal)
print("done fixing %s" % each_vcf)


N_gatk_threshold = 50
qual_gatk_threshold = 150

sample_map_qualities={}
sample_dict={}
pos_call_dict={}
dict_qual = {}
dict_map = {}
list_pass = []
list_amb = []
returned_qual = []
returned_map = []
found_positions = {}

list_of_files = glob.glob('*vcf')
for file_name in list_of_files:
    vcf_reader = vcf.Reader(open(file_name, 'r'))
    for record in vcf_reader:
        record_position = str(record.CHROM) + "-" + str(record.POS)
        absolute_positon = record_position
        if str(record.ALT[0]) != "None" and str(record.INFO['MQ']) != "nan":
            sample_map_qualities.update({record_position:record.INFO['MQ']})
        if str(record.ALT[0]) != "None" and len(record.ALT[0]) == 1 and record.INFO['AC'][0] == 2 and record.QUAL > N_gatk_threshold:
            sample_dict.update({record_position:record.ALT[0]})
        if record.ALT[0]:
            pos_call_dict.update({record.POS:record.ALT[0]})
        else:
            pos_call_dict.update({record.POS:record.REF})
        chrom = record.CHROM
        position = record.POS
        if str(record.ALT[0]) != "None" and record.INFO['AC'][0] == 2 and len(record.REF) == 1 and record.QUAL > 150:
            pass
        
        try:
            if int(record.QUAL) > 0:
                returned_qual.append(record.QUAL)
                returned_map.append(record.INFO['MQ'])
                dict_qual[absolute_positon] = returned_qual
                dict_map[absolute_positon] = returned_map
        except Exception:
            pass

        try:
            record_alt_length = len(record.ALT[0])
        except TypeError:
            record_alt_length = 0
        try:
            record_ref_length = len(record.REF)
        except TypeError:
            record_alt_length = 0

        try:
            if str(record.ALT[0]) != "None" and record_ref_length == 1 and record_alt_length == 1 and record.INFO['AC'][0] == 2 and record.QUAL > qual_gatk_threshold and record.INFO['MQ'] > 45:
                list_pass.append(absolute_positon)
            # capture ambigous defining SNPs in htmlfile
            elif str(record.ALT[0]) != "None" and record.INFO['AC'][0] == 1:
                list_amb.append(absolute_positon)
        except ZeroDivisionError:
            print ("bad line in %s at %s" % (each_vcf, absolute_positon))

########

        try:
            chrom = record.CHROM
            position = record.POS
            absolute_positon = str(chrom) + "-" + str(position)
            filter=record.FILTER
            
            # Usable positins are those that:

            # ADD PARAMETERS HERE TO CHANGE WHAT'S SNP WILL BE USED
            # IF NOT FOUND HERE THE SNP WILL BE IGNORED.  WILL NOT BE REPRESENTED.  HARD REMOVAL
            
            ## GATK parameters
            # str(record.ALT[0]) != "None" --> filter deletions
            # len(record.REF) == 1 --> filter bad ref call with 2 nt present
            # len(record.ALT[0]) == 1 --> filter bad alt call with 2 nt present
            # record.heterozygosity == 0.0 --> filter AC=1, heterozygosity.
            # record.QUAL > 150 --> filter poor quality
            # record.INFO['MQ'] --> filter low map quality
            try:
                if str(record.ALT[0]) != "None" and record.INFO['AC'][0] == 2 and len(record.REF) == 1 and record.QUAL > qual_gatk_threshold:
                    found_positions.update({absolute_positon:record.REF})
            except KeyError:
                pass
        except ZeroDivisionError:
            print ("ZeroDivisionError error found")
        except ValueError:
            print ("ValueError error found")
        except UnboundLocalError:
            print ("UnboundLocalError error found")
        except TypeError:
            print ("TypeError error found")

    print("%s passed" % file_name)