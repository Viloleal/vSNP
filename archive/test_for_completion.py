#!/usr/bin/env python

import glob
import os


#set working directory to script1 folder

write_out = open("completion_test.txt", 'w')

directory_list = glob.glob('./*')

for each_dir in directory_list:
    if glob.glob('./' + each_dir + '/alignment/*_zc.vcf'):
        print("{} is good" .format(each_dir))
        print("{} is good" .format(each_dir), file=write_out)
    else:
        print("######## \t{}\t needs to be fixed". format(each_dir))
        print("######## \t{}\t needs to be fixed". format(each_dir), file=write_out)
        
write_out.close()