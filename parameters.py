def parameters(give_option):

    real_path = os.path.dirname(os.path.realpath(__file__))
    print("real path command --> {}".format(real_path))
    real_path = real_path.split('/')
    root_path = '/'.join(real_path[:-1])
    dependents_dir = root_path + "/dependencies"
    
    if os.path.isdir("/bioinfo11/TStuber/Results"): #check bioinfo from server
        upload_to = "/bioinfo11/TStuber/Results"
    else:
        upload_to = None

    if give_option == "salmonella":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/gen-bact/salmonella/snp_pipeline/script_dependents/script1"
        upload_to, script_dependents = script1.update_directory(dependents_dir)#***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/NC_016856-NC_016855.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/NC_016856-NC_016855HighestQualitySNPs.vcf"
        gbk_file = script_dependents + "/NC_016856-NC_016855.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found

    if give_option == "ab1":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/brucella/abortus1/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/NC_00693c.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/NC_00693cHighestQualitySNPs.vcf"
        gbk_file = script_dependents + "/NC_006932-NC_006933.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found
    
    if give_option == "ab3":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/brucella/abortus3/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/CP007682-7683c.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/CP007682-7683cHighestQualitySNPs.vcf"
        gbk_file = script_dependents + "/CP007682-CP007683.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"
        

        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found
    
    if give_option == "canis":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/brucella/canis/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/BcanisATCC23365.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/canisHighestQualitySNPs.vcf"
        gbk_file = script_dependents + "/NC_010103-NC_010104.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found

    if give_option == "ceti1":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/brucella/ceti1/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/Bceti1Cudo.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/ceti1HighestQualitySNPs.vcf"
        gbk_file = None #script_dependents + "/no.gff"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found

    if give_option == "ceti2":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/brucella/ceti2/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/Bceti2-TE10759.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/ceti2HighestQualitySNPs.vcf"
        gbk_file = script_dependents + "/NC_022905-NC_022906.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found

    if give_option == "mel1":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/brucella/melitensis-bv1/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/mel-bv1-NC003317.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/mel-bv1-NC003317-highqualitysnps.vcf"
        gbk_file = script_dependents + "/NC_003317-NC_003318.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found

    if give_option == "mel1b":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/brucella/melitensis-bv1b/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/mel-bv1b-CP018508.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/mel-bv1b-CP018508-highqualitysnps.vcf"
        gbk_file = script_dependents + "/mel-bv1b-CP018508.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found

    if give_option == "mel2":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/brucella/melitensis-bv2/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/mel-bv2-NC012441.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/mel-bv2-NC012441-highqualitysnps.vcf"
        gbk_file = script_dependents + "/NC_012441-NC_012442.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found
        
    if give_option == "mel3":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/brucella/melitensis-bv3/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/mel-bv3-NZCP007760.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/mel-bv3-NZCP007760-highqualitysnps.vcf"
        gbk_file = script_dependents + "/NZ_CP007760-NZ_CP007761.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found

    if give_option == "suis1":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/brucella/suis1/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/NC_017251-NC_017250.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/B00-0468-highqualitysnps.vcf"
        gbk_file = script_dependents + "/NC_017251-NC_017250.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found

    if give_option == "suis2":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/brucella/suis2/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/NC_010169-NC_010167.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/BsuisF7-06-1-highqualitysnps.vcf"
        gbk_file = script_dependents + "/NC_010169-NC_010167.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found

    if give_option == "suis3":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/brucella/suis3/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/NZ_CP007719-NZ_CP007718.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/highqualitysnps.vcf"
        gbk_file = None #script_dependents + "/no.gff"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found

    if give_option == "suis4":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/brucella/suis4/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/B-REF-BS4-40.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/suis4HighestQualitySNPs.vcf"
        gbk_file = None #script_dependents + "/no.gff"
        email_list = "tod.p.stuber@aphis.usda.gov"

        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found
        
    if give_option == "ovis":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/brucella/ovis/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/BovisATCC25840.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/BovisATCC25840HighestQualitySNPs.vcf"
        gbk_file = script_dependents + "/NC_009505-NC_009504.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"

        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found
        
    if give_option == "neo":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/brucella/neotomae/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/KN046827.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/ERR1845155-highqualitysnps.vcf"
        gbk_file = script_dependents + "/KN046827.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"

        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found

    if give_option == "bovis":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/mycobacterium/tbc/tbbov/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/spoligotype_db.txt"
        reference = script_dependents + "/NC_002945.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/HighestQualitySNPs.vcf"
        gbk_file = script_dependents + "/NC_002945.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found
        
    if give_option == "af":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/mycobacterium/tbc/af2122/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/spoligotype_db.txt"
        reference = script_dependents + "/NC_002945v4.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/highqualitysnps.vcf"
        gbk_file = script_dependents + "/NC_002945v4.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found

    if give_option == "h37":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/mycobacterium/tbc/h37/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/spoligotype_db.txt"
        reference = script_dependents + "/NC000962.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/15-3162-highqualitysnps.vcf"
        gbk_file = script_dependents + "/NC_000962.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found

    if give_option == "para":
        found=True
        #Remove network path at and left of "Results"
        dependents_dir="/mycobacterium/avium_complex/para_cattle-bison/script_dependents/script1"
        upload_to, remote, script_dependents = script1.update_directory(dependents_dir) #***FUNCTION CALL
        
        spoligo_db = script_dependents + "/nospoligo.txt"
        reference = script_dependents + "/NC_002944.fasta"
        print("Reference being used: %s" % reference)
        hqs = script_dependents + "/HQ-NC002944.vcf"
        gbk_file = script_dependents + "/NC_002944.gbk"
        email_list = "tod.p.stuber@aphis.usda.gov"
        
        option_list=[dependents_dir, reference, hqs, gbk_file, email_list, upload_to, remote, script_dependents, spoligo_db]
        return option_list, found


################
# SCRIPT 2
################

if args.species == "salmonella":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/gen-bact/salmonella/snp_pipeline/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) # returned upload_to, remote, local (aka: script_dependents) --> local is where working dependencies are located
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/gen-bact/salmonella/snp_pipeline/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx" # this may not be available if there is no access to f drive.  f drive record will not get cp to cut bioinfo list and then cp locally.  Can also manually put something in ~/dependencies on github.
    gbk_file = script_dependents + "/NC_016856-NC_016855.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/gen-bact/salmonella/snp_pipeline/script2"
    excelinfile = script_dependents + "/Filtered_Regions.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov"
    
elif args.species == "suis1":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/suis1/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) # returned upload_to, remote, local (aka: script_dependents) --> local is where working dependencies are located
    
    bruc_private_codes(upload_to) # if f drive then upload fixed column 32 to bioinfo
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx" # this may not be available if there is no access to f drive.  f drive record will not get cp to cut bioinfo list and then cp locally.  Can also manually put something in ~/dependencies on github.
    gbk_file = script_dependents + "/NC_017251-NC_017250.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/suis1/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, christine.r.quance@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "suis2":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/suis2/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) # returned upload_to, remote, local (aka: script_dependents) --> local is where working dependencies are located
    
    bruc_private_codes(upload_to) # if f drive then upload fixed column 32 to bioinfo
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx" # this may not be available if there is no access to f drive.  f drive record will not get cp to cut bioinfo list and then cp locally.  Can also manually put something in ~/dependencies on github.
    gbk_file = script_dependents + "/NC_010169-NC_010167.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/Defining_SNPs.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/suis2/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_Suis2.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, christine.r.quance@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "suis3":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/suis3/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) # returned upload_to, remote, local  --> local is where working dependencies are located
    bruc_private_codes(upload_to)
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NZ_CP007719-NZ_CP007718.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/suis3/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, christine.r.quance@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "suis4":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/suis4/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    bruc_private_codes(upload_to)
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    #gbk_file = script_dependents + ""
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/suis4/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, christine.r.quance@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "ab1":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/abortus1/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    bruc_private_codes(upload_to)
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_006932-NC_006933.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/abortus1/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, christine.r.quance@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "ab3":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/abortus3/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    bruc_private_codes(upload_to)
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/CP007682-CP007683.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/abortus3/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, christine.r.quance@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "mel1":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/melitensis-bv1/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    bruc_private_codes(upload_to)
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_003317-NC_003318.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/melitensis-bv1/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, christine.r.quance@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "mel1b":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/melitensis-bv1b/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    bruc_private_codes(upload_to)
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/mel-bv1b-CP018508.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/melitensis-bv1b/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, christine.r.quance@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "mel2":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/melitensis-bv2/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    bruc_private_codes(upload_to)
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_012441-NC_012442.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/melitensis-bv2/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, christine.r.quance@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "mel3":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/melitensis-bv3/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    bruc_private_codes(upload_to)
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NZ_CP007760-NZ_CP007761.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/melitensis-bv3/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, christine.r.quance@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "canis":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/canis/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    bruc_private_codes(upload_to)
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_010103-NC_010104.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/canis/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, christine.r.quance@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "ceti1":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/ceti1/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    bruc_private_codes(upload_to)
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    #gbk_file = script_dependents + ""
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/ceti1/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, christine.r.quance@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "ceti2":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/ceti2/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    bruc_private_codes(upload_to)
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_022905-NC_022906.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/ceti2/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, christine.r.quance@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"
        
elif args.species == "ovis":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/ovis/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    bruc_private_codes(upload_to)
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_009505-NC_009504.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/ovis/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, christine.r.quance@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"
        
elif args.species == "neo":

    qual_gatk_threshold = 300
    N_gatk_threshold = 350
    
    #Remove network path at and left of "Results"
    dependents_dir="/brucella/neotomae/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    bruc_private_codes(upload_to)
    try:
        shutil.copy(upload_to + "/brucella/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/KN046827.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/brucella/neotomae/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions.xlsx"
    print(excelinfile)
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, christine.r.quance@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "bovis":
    
    qual_gatk_threshold = 150
    N_gatk_threshold = 150
    
    #Remove network path at and left of "Results"
    dependents_dir="/mycobacterium/tbc/tbbov/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/mycobacterium/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")
    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_002945.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/mycobacterium/tbc/tbbov/script2"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    filter_files = script_dependents + "/filter_files"
    print ("filter_files %s" % filter_files)
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:
        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if not args.email:
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "af":
    
    qual_gatk_threshold = 150
    N_gatk_threshold = 150
    
    #Remove network path at and left of "Results"
    dependents_dir="/mycobacterium/tbc/af2122/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/mycobacterium/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")
    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_002945v4.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/mycobacterium/tbc/af2122/script2"
    excelinfile = script_dependents + "/Filtered_Regions.xlsx"
    filter_files = script_dependents + "/filter_files"
    print ("filter_files %s" % filter_files)
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:
        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if not args.email:
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "h37":
    
    qual_gatk_threshold = 150
    N_gatk_threshold = 150
    
    #Remove network path at and left of "Results"
    dependents_dir="/mycobacterium/tbc/h37/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/mycobacterium/genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_000962.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/mycobacterium/tbc/h37/script2"
    excelinfile = script_dependents + "/Filtered_Regions_python.xlsx"
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

elif args.species == "para":
    
    qual_gatk_threshold = 150
    N_gatk_threshold = 150
    
    #Remove network path at and left of "Results"
    dependents_dir="/mycobacterium/avium_complex/para_cattle-bison/script_dependents/script2"
    
    upload_to, remote, script_dependents = update_directory(dependents_dir) #***FUNCTION CALL
    try:
        shutil.copy(upload_to + "/mycobacterium/avium_complex/avium_genotyping_codes.xlsx", script_dependents)
    except FileNotFoundError:
        print ("will use previously used genotyping_codes.xlsx file")

    genotypingcodes = script_dependents + "/avium_genotyping_codes.xlsx"
    gbk_file = script_dependents + "/NC_002944.gbk"
    # This file tells the script how to cluster VCFs
    definingSNPs = script_dependents + "/DefiningSNPsGroupDesignations.xlsx"
    remove_from_analysis = script_dependents + "/RemoveFromAnalysis.xlsx"
    bioinfoVCF = upload_to + "/mycobacterium/avium_complex/para_cattle-bison/vcfs"
    excelinfile = script_dependents + "/Filtered_Regions.xlsx"
    filter_files = script_dependents + "/filter_files"
    if os.path.isdir(filter_files):
        shutil.rmtree(filter_files)
        os.mkdir(filter_files)
    else:        os.mkdir(filter_files)
    get_filters(excelinfile, filter_files) #***FUNCTION CALL
    if args.email == "s":
        email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov"

else:
    parser.print_help()
    print ("\n#####EXIT AT SETTING OPTIONS, Check that a \"-s\" species was provided\n")
    sys.exit(0)