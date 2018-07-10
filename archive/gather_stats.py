#!/usr/bin/env python

import glob
import xlrd

excel_stat_files = glob.glob("*xlsx")

write_out = open ("stats_collection.txt", 'w')

for excelinfile in excel_stat_files:
    wb = xlrd.open_workbook(excelinfile)
    ws = wb.sheet_by_index(0) #in first worksheet
    data_row = ws.row(1) #list
    for mydata in data_row:
        print(mydata.value, end='\t', file=write_out) #mydata is type: Cell value
    print("", file=write_out)

write_out.close()