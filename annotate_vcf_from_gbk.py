import pandas as pd


def annotate_vcf(vcf_path_in, gbk_path_in):
    print("vcf_path_in: {}".format(vcf_path_in))
    print("gbk_path_in: {}".format(gbk_path_in))
    df = pd.read_csv(vcf_path_in, sep='\t', names=["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "sample"], comment='#')
    print(df.head())
    

vcf_path_in = '/Users/tstuber/Desktop/vcf_test_files/bovis/test/00-11MIDNR_zc.vcf'
gbk_path_in = '/Users/tstuber/workspace/vSNP/dependencies/mycobacterium/tbc/af2122/script_dependents/NC_002945v4.gbk'
annotate_vcf(vcf_path_in, gbk_path_in)