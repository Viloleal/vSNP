###############################################
###############################################
##################script2######################
###############################################
###############################################

    home = os.path.expanduser("~")
    
def setup_raxml():    
    # IF AVX2 IS AVAILABE (CHECK WITH `cat /proc/cpuinfo | grep -i "avx"`). CREATE A LINK TO: `ln -s path_to_raxmlHPC-PTHREADS-AVX2 raxml.  Place "raxml" in your path.  This will allow "raxml" to be found first which will call AVX2 version of RAxML
    
    try:
        subprocess.call("raxml", stdout=open(os.devnull, 'wb'))
        sys_raxml = "raxml"
        #print ("%s found" % sys_raxml)
    except OSError:
        print ("looking for RAxML")
        try:
            subprocess.call("raxmlHPC-PTHREADS")
            sys_raxml = "raxmlHPC-PTHREADS"
            print ("%s found" % sys_raxml)
        except OSError:
            try:
                subprocess.call("raxmlHPC-SSE3")
                sys_raxml = "raxmlHPC-SSE3"
                print ("%s found" % sys_raxml)
            except OSError:
                print ("looking for RAxML")
                try:
                    subprocess.call("raxmlHPC")
                    sys_raxml = "raxmlHPC"
                    print ("RAxML found")
                except OSError:
                    print ("#####RAxML is not in you PATH")
                    print ("#####See help page for support")
                    sys.exit(0)

    print ("\n\n----> RAxML found in $PATH as: %s <-----" % sys_raxml)
    
def update_directory(dependents_dir): # UPDATE DIRECTORIES
    home = os.path.expanduser("~")
    print("dependents_dir %s\n" % dependents_dir)
    
    if os.path.isdir("/bioinfo11/TStuber/Results"): #check bioinfo from server
        upload_to = "/bioinfo11/TStuber/Results"
        remote="/bioinfo11/TStuber/Results" + dependents_dir
        if os.path.isdir("/Users/Shared"):
            dep_path = "/Users/Shared"
            dir_split = dependents_dir.split('/')[1:]
            for i in dir_split:
                dep_path += '/' + i
                if not os.path.exists(dep_path):
                    os.makedirs(dep_path)
            local = "/Users/Shared" + dependents_dir
            if os.path.isdir(local):
                try:
                    shutil.rmtree(local)
                    shutil.copytree(remote, local)
                except:
                    pass
        elif os.path.isdir("/home/shared"):
            dep_path = "/home/shared"
            dir_split = dependents_dir.split('/')[1:]
            for i in dir_split:
                dep_path += '/' + i
                if not os.path.exists(dep_path):
                    os.makedirs(dep_path)
            local = "/home/shared" + dependents_dir
            if os.path.isdir(local):
                try:
                    shutil.rmtree(local)
                    shutil.copytree(remote, local)
                except:
                    pass

    elif os.path.isdir("/Volumes/root/TStuber/Results"): #check bioinfo from Mac
        upload_to = "/Volumes/root/TStuber/Results"
        remote="/Volumes/root/TStuber/Results" + dependents_dir
        if os.path.isdir("/Users/Shared"):
            dep_path = "/Users/Shared"
            dir_split = dependents_dir.split('/')[1:]
            for i in dir_split:
                dep_path += '/' + i
                if not os.path.exists(dep_path):
                    os.makedirs(dep_path)
            local = "/Users/Shared" + dependents_dir
            if os.path.isdir(local):
                try:
                    shutil.rmtree(local)
                    shutil.copytree(remote, local)
                except:
                    pass
        elif os.path.isdir("/home/shared"):
            dep_path = "/home/shared"
            dir_split = dependents_dir.split('/')[1:]
            for i in dir_split:
                dep_path += '/' + i
                if not os.path.exists(dep_path):
                    os.makedirs(dep_path)
            local = "/home/shared" + dependents_dir
            if os.path.isdir(local):
                try:
                    shutil.rmtree(local)
                    shutil.copytree(remote, local)
                except:
                    pass

    #### PLACE A CHECK FROM GITHUB
    
    elif os.path.isdir("/Users/Shared" + dependents_dir): #check local copy in shared folder
        upload_to ="not_found"
        remote = "not_found"
        local = "/Users/Shared" + dependents_dir
            
    elif os.path.isdir("/home/shared" + dependents_dir): #check local copy in shared folder
        upload_to ="not_found"
        remote = "not_found"
        local = "/home/shared" + dependents_dir
    
    elif os.path.isdir(home + "/dependencies" + dependents_dir): #check local copy from Git repo
        upload_to ="not_found"
        remote = "no remote"
        script_location = home # points to home directory
        local = home + "/dependencies" + dependents_dir # sets dependencies directory to home directory
    else:
        os.makedirs(home + "/dependencies")
        print("\n\nDOWNLOADING DEPENDENCIES FROM GITHUB... ***\n\n")
        git.Repo.clone_from("https://github.com/USDA-VS/dependencies.git", home + "/dependencies")
        upload_to ="not_found"
        remote = "no remote"
        script_location = home # points to home directory
        local = home + "/dependencies" + dependents_dir # sets dependencies directory to home directory
    
    print("\n####################DIRECTORY LOCATION")
    print("####################upload_to: %s" % upload_to)
    print("####################remote: %s" % remote)
    print("####################local: %s\n" % local)
    
    return upload_to, remote, local

# Get filters set up
def get_filters(excelinfile, filter_files):
    for i in glob.glob(filter_files + "/*"):
        os.remove(i)

    wb = xlrd.open_workbook(excelinfile)
    sheets = wb.sheet_names()
    for sheet in sheets:
        ws = wb.sheet_by_name(sheet)

        myrange = lambda start, end: range(start, end+1)

        for colnum in range(ws.ncols): # for each column in worksheet
            file_out = filter_files + "/" + ws.col_values(colnum)[0] + ".txt" # column header naming file
            write_out = open (file_out, 'at')
            mylist = ws.col_values(colnum)[1:] # list of each field in column, minus the header
            mylist = [x for x in mylist if x] # remove blank cells
            for value in mylist:
                value = str(value)
                value = value.replace(sheet + "-", '')
                if "-" not in value:
                    value=int(float(value)) # change str to float to int
                    print (sheet + "-" + str(value), file=write_out)
                elif "-" in value:
                    value = value.split("-")
                    for i in range(int(value[0]), int(value[1]) + 1 ):
                        print (sheet + "-" + str(i), file=write_out)
    write_out.close()

##################
# FUNCTIONS
##################

# Test for duplicate samples
def test_duplicate():
    dup_list = []
    list_of_files = glob.glob('*vcf')
    for line in list_of_files:
        line=re.sub(r'(.*)[_.].*', r'\1', line)
        dup_list.append(line)
    # find duplicates in list
    duplicates = [k for k,v in Counter(dup_list).items() if v>1]
    if len(duplicates) > 0:
        print ("Duplicates Found: %s " % duplicates)
        print ("\n***Error:  Duplicate VCFs")
        sys.exit(0)
    else:
        print ("\nno duplicate VCFs\n")

# Change file names
def change_names():
    global malformed
    code_dictionary = {}
    try:
        wb = xlrd.open_workbook(genotypingcodes)
        ws = wb.sheet_by_index(0)
        for rownum in range(ws.nrows):
            new_name = str(ws.row_values(rownum)[0])
            new_name = new_name.rstrip()
            new_name = re.sub('[\/() ]', '_', new_name)
            new_name = re.sub('#', 'num', new_name)
            new_name = re.sub('_-', '_', new_name)
            new_name = re.sub('-_', '_', new_name)
            new_name = re.sub('__+', '_', new_name)
            new_name = re.sub('_$', '', new_name)
            new_name = re.sub('-$', '', new_name)
            new_name = re.sub(',', '', new_name)
            try:
                elite_test = ws.row_values(rownum)[1]
            except IndexError:
                #print ("except IndexError: when changing names")
                elite_test = ""
            #print("newname %s" % new_name)
            try:
                if new_name[-1] != "_":
                    new_name = new_name + "_"
            except IndexError:
                pass
            code_dictionary.update({new_name:elite_test})
    except FileNotFoundError:
        print ("\n#### except: FileNotFoundError, there was not a \"genotypingcodes\" file given to change names\n")

    names_not_changed = []
    list_of_files = glob.glob('*vcf')
    for each_vcf in list_of_files:
        vcf_found = False
        vcf_pretext = re.sub(r'(.*?)[._].*', r'\1', each_vcf) # ? was needed to make greedy, in my view the regex was searching right to left without it.
        vcf_pretext = vcf_pretext.rstrip()
        #Added '^' because h37 18-2397 was finding bovis 18-011018-2397, 2018-06-19
        myregex = re.compile('^' + vcf_pretext + '_.*') #underscore required to make myregex.search below greedy.  so it finds exact match and not all matches. ex: 10-01 must match 10-01 not 10-010 also
        for k, v in code_dictionary.items():
            try:
                if myregex.search(k):
                    k= k.strip('_')
                    #print("myregex %s, matches %s" % (myregex, k))
                    os.rename(each_vcf, k + ".vcf")
                    vcf_found = True
            except FileNotFoundError:
                print ("except FileNotFoundError %s" % each_vcf)
        if vcf_found == False:
                    names_not_changed.append(each_vcf)
    names_not_changed = set(names_not_changed) # remove duplicates

    if args.elite:
        list_of_files = []
        list_of_files = glob.glob('*vcf')
        if not os.path.exists("temp_hold"):
            print ("making temp_hold directory")
            os.makedirs("temp_hold") # make all_vcfs if none exists
        for each_vcf in list_of_files:
            # Default 1 * 24 * 60 *60
            time_test = time.time() - os.path.getmtime(each_vcf) < (1 * 24 * 60 *60) # 1day * (24*60*60)sec in day
            print ("%s each_vcf" % each_vcf)
            vcf_pretext = re.sub(r'(.*?)[._].*', r'\1', each_vcf) # ? was needed to make greedy, in my view the regex was searching right to left without it.
            vcf_pretext = vcf_pretext.rstrip()
            myregex = re.compile(vcf_pretext + '.*')
            if time_test:
                print ("time_test true %s" % each_vcf)
                shutil.copy(each_vcf, "temp_hold")
            else:
                for k, v in code_dictionary.items():
                    if myregex.search(k):
                        try:
                            print ("##### %s" % time_test)
                            if v == "Yes": # if marked yes in column 2 of genotyping codes
                                print ("marked yes %s" % each_vcf)
                                shutil.copy(each_vcf, "temp_hold") # if "Yes" then moved to temp_hold
                            else:
                                print ("file will be discarded %s" % each_vcf)
                        except FileNotFoundError:
                            print ("except FileNotFoundError %s" % each_vcf)
                os.remove(each_vcf)
        shutil.rmtree('starting_files')
        os.makedirs('starting_files')
        os.renames('temp_hold', 'starting_files')
        list_of_files = glob.glob('starting_files/*vcf')
        file_number = len(list_of_files) # update the file_number to present on summary
        for each_vcf in list_of_files:
            shutil.copy(each_vcf, root_dir)
        all_starting_files = glob.glob('*vcf')
        print (file_number)
    
    #fix files
    vcf_list = glob.glob('*vcf')
    print("Fixing files...\n")
    if args.debug_call and not args.get:
        for each_vcf in vcf_list:
            print(each_vcf)
            mal = fix_vcf(each_vcf)
            malformed = malformed + list(mal)
    else:
        with Pool(maxtasksperchild=4) as pool:
            mal = pool.map(fix_vcf, vcf_list, chunksize=8)
            malformed = malformed + list(mal)
    print("done fixing")

    return names_not_changed

def send_email():
    print ("Sending Email...")
    print ("Sending to:")

    msg = MIMEMultipart()
    msg['From'] = "tod.p.stuber@aphis.usda.gov"
    msg['To'] = email_list
    msg['Subject'] = "Script 2 " + args.species
    with open(htmlfile_name) as fp:
        msg.attach(MIMEText(fp.read(), 'html'))

    part = MIMEBase('application', "octet-stream")
    part.set_payload(open("summary_log.html", "r").read())
    encoders.encode_base64(part)
    part.add_header('Content-Disposition', 'attachment; filename="summary_log.html"')
    msg.attach(part)

    smtp = smtplib.SMTP('10.10.8.12')
    smtp.send_message(msg)

    #smtp.send_message(msg)
    #smtp.send_message(msg.as_string())
    #smtp.sendmail(email_list, msg.as_string())
    #smtp.sendmail("tod.p.stuber@aphis.usda.gov", email_list, msg.as_string())
    smtp.quit()

    if args.email == "none":
        print ("\n\temail not sent")
    elif args.email:
        send_email()
        print ("\n\temail sent to: %s" % email_list)
    else:
        print ("\n\temail not sent")

    if args.upload:
        print ("Uploading Samples...")
        def copytree(src, dst, symlinks=False, ignore=None): #required to ignore permissions
            try:
                for item in os.listdir(src):
                    s = os.path.join(src, item)
                    d = os.path.join(dst, item)
                    try:
                        if os.path.isdir(s):
                            shutil.copytree(s, d, symlinks, ignore)
                        else:
                            shutil.copy2(s, d)
                    except shutil.Error:
                        pass
            except FileNotFoundError:
                print ("except FileNotFoundError: file not found")

        #upload to bioinfoVCF
        src = root_dir
        dst = bioinfoVCF + "/" + os.path.basename(os.path.normpath(root_dir))
        print ("\n\t%s is copying to %s" % (src, dst))
        os.makedirs(dst, exist_ok=True)
        copy_tree(src, dst, preserve_mode=0, preserve_times=0)

def flatten(l):
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el

def get_pretext_list(in_list):
    outlist = []
    for i in in_list:
        pretext  = re.sub('[_.].*', '', i)
        outlist.append(pretext)
    return outlist

def get_snps(directory):
    os.chdir(root_dir+ "/" + directory)
    print ("\n----------------------------")
    print ("\nworking on: %s " % directory)
    outdir=str(os.getcwd()) + "/"
    # FILTER position all list
    list_filter_files = glob.glob(filter_files + '/*')

    filter_file = "empty" # if filter an all_vcf file not found mark as empty
    filter_group = "empty" # if a group specific filter file is not found mark as empty
    for i in list_filter_files:
        if "-All.txt" in i:
            filter_file = i

    for i in list_filter_files:
        if directory  + ".txt" in i:
            filter_group = i

    print ("%s --> filter_file %s " % (directory, filter_file))
    print ("%s --> filter_group %s " % (directory, filter_group))
    print ("%s --> outdir %s " % (directory, outdir))

    files = glob.glob('*vcf')
    all_positions = {}
    if args.debug_call and not args.get:
        for i in files:
            found_positions = find_positions(i)
            all_positions.update(found_positions)
    else:
        with Pool(maxtasksperchild=4) as pool:
            for found_positions in pool.map(find_positions, files, chunksize=8):
                all_positions.update(found_positions)

    print ("Directory %s found positions %s" % (directory, len(all_positions)))
    presize=len(all_positions)

    # Filter applied to all positions
    if not filter_file is "empty":
        with open(filter_file, 'rt') as f:
            filter_list = f.read().splitlines() #removes \n
        for pos in filter_list:
            all_positions.pop(pos, None)
        f.close()

    # Filter applied to group
    if not filter_group is "empty":
        with open(filter_group, 'rt') as f:
            filter_list = f.read().splitlines() #removes \n
        for pos in filter_list:
            all_positions.pop(pos, None)
        f.close()

    print ("\nDirectory: ", directory)
    print ("Total positions found: %s" % format(presize, ",d"))
    print ("Possible positions filtered %s" % format(len(filter_list), ",d"))
    print ("Positions after filtering %s\n" % format(len(all_positions), ",d"))

    if args.filter:
        #write to files
        positions_to_filter = "positions_to_filter.txt"
        positions_to_filter_details = "positions_to_filter_details.txt"
        good_snps = "good_snps_details.txt"
        write_out_positions=open(positions_to_filter, 'w')
        write_out_details=open(positions_to_filter_details, 'w')
        write_out_good_snps=open(good_snps, 'w')

        files = glob.glob('*vcf')

        #calculate mean/max qual and map at all possible positions
        dd_qual = {}
        dd_map = {}
        if args.debug_call:
            for each_vcf in files:
                print ("working on: %s" % each_vcf)
                dict_qual, dict_map = find_filter_dict(each_vcf)
                keys = set(dd_qual).union(dict_qual)
                no = []
                #make position (key) and qual/maps list (value)
                dd_qual = dict((k, dd_qual.get(k, no) + dict_qual.get(k, no)) for k in keys)
                keys = set(dd_map).union(dict_map)
                no = []
                dd_map = dict((k, dd_map.get(k, no) + dict_map.get(k, no)) for k in keys)
        else:
            with Pool(maxtasksperchild=4) as pool:
                for dict_qual, dict_map in pool.map(find_filter_dict, files, chunksize=8):
                    keys = set(dd_qual).union(dict_qual)
                    no = []
                    dd_qual = dict((k, dd_qual.get(k, no) + dict_qual.get(k, no)) for k in keys)
                    keys = set(dd_map).union(dict_map)
                    no = []
                    dd_map = dict((k, dd_map.get(k, no) + dict_map.get(k, no)) for k in keys)

        #dict_qual=dict((k, v) for k, v in dict_qual.items() if v)
        #dict_map=dict((k, v) for k, v in dict_map.items() if v)

        ave_qual = {}
        max_qual = {}
        for k, v in dd_qual.items():
            #only use if > 3 positions have been called
            if len(v) > 3:
                ave_qual[k]=np.mean(v)
                max_qual[k]=np.max(v)

        #provides dictionary as key -> absolute poisiton, value -> average qual/map
        ave_map = {}
        max_map = {}
        for k, v in dd_map.items():
            if len(v) > 3:
                ave_map[k]=np.mean(v)
                max_map[k]=np.max(v)		

        # get all possible used positions
        all_maybe_filter = []
        for k in ave_qual.keys():
            all_maybe_filter.append(k)
        for k in max_qual.keys():
            all_maybe_filter.append(k)
        for k in ave_map.keys():
            all_maybe_filter.append(k)
        for k in max_map.keys():
            all_maybe_filter.append(k)
            # remove duplicates
            all_maybe_filter = list(set(all_maybe_filter))

        #remove those in filter list
        #Filter applied to all positions
        if not filter_file is "empty":
            with open(filter_file, 'rt') as f:
                filter_list = f.read().splitlines() #removes \n
                try:
                    for pos in filter_list:
                        all_maybe_filter.pop(pos)
                except TypeError:
                    pass
                except KeyError:
                    pass
            f.close()

        # Filter applied to group
        if not filter_group is "empty":
            with open(filter_group, 'rt') as f:
                filter_list = f.read().splitlines() #removes \n
                try:
                    for pos in filter_list:
                        all_maybe_filter.pop(pos)
                except TypeError:
                    pass
                except KeyError:
                    pass
            f.close()

        # for each possible posible position check if to filter.
        for absolute_positon in all_maybe_filter:
            ave_qual_value = ave_qual[absolute_positon]
            max_qual_value = max_qual[absolute_positon]
            ave_map_value = ave_map[absolute_positon]
            max_map_value = max_map[absolute_positon]
            print ("%s, max_qual_value: %s, ave_qual_value: %s, max_map_value: %s, ave_map_value: %s" % (absolute_positon, max_qual_value, ave_qual_value, max_map_value, ave_map_value))
            if max_qual_value < 1300 and ave_qual_value < 800 or ave_map_value < 56:
                print ("%s, max_qual_value: %s, ave_qual_value: %s, max_map_value: %s, ave_map_value: %s" % (absolute_positon, max_qual_value, ave_qual_value, max_map_value, ave_map_value), file=write_out_details)
                print (absolute_positon, file=write_out_positions)
            else:
                print ("%s, max_qual_value: %s, ave_qual_value: %s, max_map_value: %s, ave_map_value: %s" % (absolute_positon, max_qual_value, ave_qual_value, max_map_value, ave_map_value), file=write_out_good_snps)
        write_out_positions.close()
        write_out_details.close()
        write_out_good_snps.close()

def get_annotations_table(parsimony_positions):
    print ("Getting annotations...")
    dict_annotation = {}
    gbk_dict = SeqIO.to_dict(SeqIO.parse(gbk_file, "genbank"))
    for each_absolute_pos in parsimony_positions:
        each_absolute_pos = each_absolute_pos.split("-")
        chrom = each_absolute_pos[0]
        pos = int(each_absolute_pos[1])
        pos_found = False
        for each_key, each_value in gbk_dict.items():
            if chrom == each_key: # need to check chrom when multiple chroms present
                for feature in each_value.features:
                    if pos in feature and "CDS" in feature.type:
                        myproduct = "none list"
                        mylocus = "none list"
                        mygene = "none list"
                        myproduct = feature.qualifiers['product'][0]
                        mylocus = feature.qualifiers['locus_tag'][0]
                        if "gene" in feature.qualifiers:
                            mygene = feature.qualifiers['gene'][0]
                        myout = myproduct + ", gene: " + mygene + ", locus_tag: " + mylocus
                        pos_found = True
            if pos_found == False:
                myout = "No annotated product"
            dict_annotation.update({chrom + "-" + str(pos):myout})
    return (dict_annotation)

def get_read_mean(rec):
    mean_q = int(mean(rec.letter_annotations['phred_quality']))
    return mean_q

#map pooled from script 1
def read_aligner(directory):
    os.chdir(directory)
    R1 = glob.glob('*_R1*fastq.gz')
    R2 = glob.glob('*_R2*fastq.gz')

    ###
    read_quality_stats = {}
    print("Getting mean for {}" .format(R1[0]))
    handle = gzip.open(R1[0], "rt")
    mean_quality_list=[]
    for rec in SeqIO.parse(handle, "fastq"):
        mean_q = get_read_mean(rec)
        mean_quality_list.append(mean_q)

    read_quality_stats["Q_ave_R1"] = "{:.1f}" .format(mean(mean_quality_list))
    thirty_or_greater_count = sum(i > 29 for i in mean_quality_list)
    read_quality_stats["Q30_R1"] = "{:.1%}" .format(thirty_or_greater_count/len(mean_quality_list))

    print("Getting mean for {}" .format(R2[0]))
    handle = gzip.open(R2[0], "rt")
    mean_quality_list=[]
    for rec in SeqIO.parse(handle, "fastq"):
        mean_q = get_read_mean(rec)
        mean_quality_list.append(mean_q)

    read_quality_stats["Q_ave_R2"] = "{:.1f}" .format(mean(mean_quality_list))
    thirty_or_greater_count = sum(i > 29 for i in mean_quality_list)
    read_quality_stats["Q30_R2"] = "{:.1%}" .format(thirty_or_greater_count/len(mean_quality_list))
    ###

    if args.species:
        sample = script1(R1[0], R2[0], args.species) #force species
    else:
        sample = script1(R1[0], R2[0]) #no species give, will find best
    try:
        stat_summary = sample.align_reads(read_quality_stats)
        return(stat_summary)
        for k, v in stat_summary.items():
            print("%s: %s" % (k, v))
    except:
        print("### Unable to return stat_summary")
        return #(stat_summary)
        pass

def fix_vcf(each_vcf):
    mal = []
    ###
    # Fix common VCF errors
    if args.debug_call and not args.get:
        print ("FIXING FILE: " + each_vcf)
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

def find_filter_dict(each_vcf):
    dict_qual = {}
    dict_map = {}
    vcf_reader = vcf.Reader(open(each_vcf, 'r'))
    for record in vcf_reader:
        absolute_positon = str(record.CHROM) + "-" + str(record.POS)
        try:
            returned_qual = []
            returned_map = []
            if int(record.QUAL) > 0:
                returned_qual.append(record.QUAL)
                returned_map.append(record.INFO['MQ'])
                dict_qual[absolute_positon] = returned_qual
                dict_map[absolute_positon] = returned_map
        except Exception:
            pass
    return dict_qual, dict_map
















print ("\nSET VARIABLES")
print ("\tgenotypingcodes: %s " % genotypingcodes)







htmlfile_name = root_dir+ "/summary_log.html"
htmlfile = open(htmlfile_name, 'at')

startTime = datetime.now()
print ("\n\n*** START ***\n")
print ("Start time: %s" % startTime)

# DIRECTORY TEST AND BACKUP
if getattr(sys, 'frozen', False):
    script_used = os.path.realpath(sys.executable)
elif __file__:
    script_used = os.path.realpath(__file__)

print ("\nScript used: %s \n" % script_used)

# make backup
os.makedirs('starting_files')
all_starting_files = glob.glob('*vcf')
for i in all_starting_files:
    shutil.copy(i, 'starting_files')



test_duplicate() #***FUNCTION CALL

global mygbk
try:
    mygbk = True
    print ("\tgbk_file: %s " % gbk_file)
except NameError:
    mygbk = False
    print ("There is not a gbk file available")
print ("\tdefiningSNPs: %s " % definingSNPs)
print ("\texcelinfile: %s " % excelinfile)
print ("\tremove_from_analysis: %s " % remove_from_analysis)
print ("\tfilter_files: %s " % filter_files)
print ("\tbioinfoVCF: %s \n" % bioinfoVCF)
###

if os.path.isfile(genotypingcodes):
    print ("\nChanging the VCF names")
    names_not_changed = change_names() # check if genotypingcodes exist.  if not skip.
else:
    print("No mapping file for VCF names")
    names_not_changed = glob.glob("*.vcf")

files = glob.glob('*vcf')
print ("REMOVING FROM ANALYSIS...")
wb = xlrd.open_workbook(remove_from_analysis)
ws = wb.sheet_by_index(0)
for each_sample in ws.col_values(0):
    each_sample = str(each_sample)
    each_sample = re.sub(r'(.*?)[._].*', r'\1', each_sample)
    #print("each sample %s" % each_sample)
    myregex = re.compile(each_sample + '.*') # create regular expression to search for in VCF list
    #print("myregex %s" % myregex)
    for i in files:
        if myregex.search(i):
            print ("### --> %s removed from the analysis" % i)
            #print (files)
            #print ("\n<h4>### --> %s removed from the analysis</h4>" % i, file=htmlfile)
            try:
                os.remove(i)
            except FileNotFoundError:
                print ("FileNotFoundError:")
vcf_starting_list = glob.glob("*.vcf")

print ("CHECKING FOR EMPTY FILES...")
files = glob.glob('*vcf')
for i in files:
    if os.stat(i).st_size == 0:
        print ("### %s is an empty file and has been deleted" % i)
        malformed.append("File was empty %s" % i)
        os.remove(i)

all_starting_files = glob.glob('*vcf')
file_number = len(all_starting_files)

print ("SORTING FILES...")
global defining_snps
defining_snps = {}
global inverted_position
inverted_position = {}
wb = xlrd.open_workbook(definingSNPs)
ws = wb.sheet_by_index(0)

for rownum in range(ws.nrows):
    position = ws.row_values(rownum)[1:][0]
    grouping = ws.row_values(rownum)[:1][0]
    # inverted positions will NOT be found in the passing positions
    # inverted positions are indicated in Defining SNPs by ending with "!"
    if position.endswith('!'):
        position = re.sub('!', '', position)
        inverted_position.update({position:grouping})
    else:
        defining_snps.update({position:grouping})
files = glob.glob('*vcf')

all_list_amb = {}
group_calls_list = []

print ("Grouping files...")
if args.debug_call and not args.get:
    for i in files:
        dict_amb, group_calls, mal = group_files(i)
        all_list_amb.update(dict_amb)
        group_calls_list.append(group_calls)
        malformed.append(mal)
else:
    with Pool(maxtasksperchild=4) as pool:
        for dict_amb, group_calls, mal in pool.map(group_files, files, chunksize=8):
            all_list_amb.update(dict_amb)
            group_calls_list.append(group_calls) # make list of list
            malformed.append(mal)
malformed = [x for x in malformed if x] # remove empty sets from list

print ("Getting directory list\n")
directory_list = next(os.walk('.'))[1] # get list of subdirectories
directory_list.remove('starting_files')

samples_in_output = []
print ("Getting SNPs in each directory")
if args.debug_call:
    for i in directory_list:
        samples_in_fasta = get_snps(i)
        samples_in_output.append(samples_in_fasta)
else:
    with futures.ProcessPoolExecutor() as pool:
        for samples_in_fasta in pool.map(get_snps, directory_list, chunksize=8):
            samples_in_output.append(samples_in_fasta)

flattened_list = []
for i in flatten(samples_in_output):
    flattened_list.append(i)
flattened_list = set(flattened_list)

count_flattened_list = len(flattened_list)
count_vcf_starting_list = len(vcf_starting_list)
start_end_file_diff_count = count_vcf_starting_list - count_flattened_list

pretext_flattened_list = get_pretext_list(flattened_list)
pretext_vcf_starting_list = get_pretext_list(vcf_starting_list)
pretext_vcf_starting_list = set(pretext_vcf_starting_list)
pretext_flattened_list.remove('root')
difference_start_end_file = pretext_vcf_starting_list.symmetric_difference(pretext_flattened_list)
difference_start_end_file = list(difference_start_end_file)
difference_start_end_file.sort()

# Zip dependency files
dependents_dir = root_dir + "/dependents"
os.makedirs(dependents_dir)
shutil.copy(definingSNPs, dependents_dir)
shutil.copy(excelinfile, dependents_dir)
zip(dependents_dir, dependents_dir)
shutil.rmtree(dependents_dir)

runtime = (datetime.now() - startTime)
print ("\n\nruntime: %s:  \n" % runtime)

#############################################
#MAKE HTML FILE:
print ("<html>\n<head><style> table { font-family: arial, sans-serif; border-collapse: collapse; width: 40%; } td, th { border: 1px solid #dddddd; padding: 4px; text-align: left; font-size: 11px; } </style></head>\n<body style=\"font-size:12px;\">", file=htmlfile)
print ("<h2>Script ran using <u>%s</u> variables</h2>" % args.species.upper(), file=htmlfile)
print ("<h4>There are %s VCFs in this run</h4>" % file_number, file=htmlfile)

#OPTIONS
print ("Additional options ran: email: %s, args.filter: %s, all_vcf: %s, elite: %s, no annotation: %s, debug: %s, get: %s, uploaded: %s" % (args.email, args.filter, args.all_vcf, args.elite, args.no_annotation, args.debug_call, args.get, args.upload), file=htmlfile)
if args.all_vcf:
    print ("\n<h4>All_VCFs is available</h4>", file=htmlfile)
elif args.elite:
    print ("\n<h4>Elite VCF comparison available</h4>", file=htmlfile)

#TIME
print ("\n<h4>Start time: %s <br>" % startTime, file=htmlfile)
print ("End time: %s <br>" % datetime.now(), file=htmlfile)
print ("Total run time: %s: </h4>" % runtime, file=htmlfile)

# ERROR LIST
if len(malformed) < 1:
    print ("<h2>No corrupt VCF removed</h2>", file=htmlfile)

else:
    print ("\n<h2>Corrupt VCF removed</h2>", file=htmlfile)
    for i in malformed:
        print ("%s <br>" % i, file=htmlfile)
    print ("<br>", file=htmlfile)

# AMBIGIOUS DEFINING SNPS
if len(all_list_amb) < 1:
    print ("\n<h2>No ambiguous defining SNPs</h2>", file=htmlfile)
else:
    print ("\n<h2>Defining SNPs are ambiguous.  They may be mixed isolates.</h2>", file=htmlfile)
    print ("<table>", file=htmlfile)
    print ("<tr align=\"left\"><th>Sample Name</th><th>Division</th><th>Absolute Position</th><tr>", file=htmlfile)
    ordered_all_list_amb = OrderedDict(sorted(all_list_amb.items()))
    for k, v in ordered_all_list_amb.items():
        k_split = k.split('\t')
        print ("<tr><td>%s</td><td>%s</td><td>%s</td></tr>" % (k_split[0], k_split[1], v), file=htmlfile)
    print ("</table>", file=htmlfile)
    print ("<br>", file=htmlfile)

#GROUPING TABLE
print ("<h2>Groupings</h2>", file=htmlfile)
print ("<table>", file=htmlfile)
print ("<tr align=\"left\"><th>Sample Name</th><tr>", file=htmlfile)

group_calls_list = list(filter(None, group_calls_list))
try:
    group_calls_list.sort(key=lambda x: x[0]) # sort list of list by first element
except IndexError:
    print("Unable to sort grouping list")
    pass

for i in group_calls_list:
    print ("<tr>", file=htmlfile)
    for x in i:
        print ("<td>%s</td>" % x, end='\t', file=htmlfile)
    print ("</tr>", file=htmlfile)
print ("</table>", file=htmlfile)

# REPORT DIFFERENCES BETWEEN STARTING FILES AND ENDING FILES REPRESENTED IN ALIGNMENTS AND TABLES
if start_end_file_diff_count < 1:
    print ("\n<h2>No files dropped from the analysis.  Input files are equal to those represented in output.</h2>", file=htmlfile)
else:
    print ("\n<h2>{} files have been dropped.  They either need a group, mixed and not finding a group or an error occured.</h2>" .format(start_end_file_diff_count), file=htmlfile)
    print ("<table>", file=htmlfile)
    print ("<tr align=\"left\"><th>Sample Name</th><tr>", file=htmlfile)
    for i in difference_start_end_file:
        print ("<tr><td>{}</td></tr>" .format(i), file=htmlfile)
    print ("</table>", file=htmlfile)
    print ("<br>", file=htmlfile)

#Capture program versions for step 2
try:
    print ("\n<h2>Program versions:</h2>", file=htmlfile)
    versions = os.popen('conda list biopython | grep -v "^#"; \
    conda list numpy | egrep -v "^#|numpydoc"; \
    conda list pandas | grep -v "^#"; \
    conda list raxml | grep -v "^#"').read()
    versions = versions.split('\n')
    for i in versions:
        print ("%s<br>" % i, file=htmlfile)
except:
    pass

#FILES NOT RENAMED
if names_not_changed:
    print ("\n<h2>File names did not get changed:</h2>", file=htmlfile)
    for i in sorted(names_not_changed):
        print ("%s<br>" % i, file=htmlfile)

print ("</body>\n</html>", file=htmlfile)
#############################################
os.chdir(root_dir)
zip("starting_files", "starting_files") # zip starting files directory
shutil.rmtree("starting_files")

htmlfile.close()




###############################################
###############################################
###################map pooled##################
###############################################
###############################################

# Group files, map pooled from script 2
def group_files(each_vcf):
    mal = ""
    list_pass = []
    list_amb = []
    dict_amb = {}
    group_calls = []
    passing = True
    #print("qual_gatk_threshold: %s " % qual_gatk_threshold)

    try:
        vcf_reader = vcf.Reader(open(each_vcf, 'r'))
        ### PUT VCF NAME INTO LIST, capturing for htmlfile
        group_calls.append(each_vcf)
            # for each single vcf getting passing position
        for record in vcf_reader:
            chrom = record.CHROM
            position = record.POS
            absolute_positon = str(chrom) + "-" + str(position)
            # find quality SNPs and put absolute positions into list
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

        for key in inverted_position.keys():
            if key not in list_pass:
                print ("key %s not in list_pass" % key)
                directory = inverted_position[key]
                print("*** INVERTED POSITION FOUND *** PASSING POSITION FOUND: \t%s\t\t%s" % (each_vcf, directory))
                if not os.path.exists(directory):
                    try:
                        os.makedirs(directory)
                    except FileExistsError:
                        null = "null"
                shutil.copy(each_vcf, directory)
                ### ADD GROUP TO LIST
                group_calls.append(directory)

        #if passing:
        # if a passing position is in the defining SNPs
        for passing_position in list_pass:
            # normal grouping
            if passing_position in defining_snps:
                directory = defining_snps[passing_position]
                print("PASSING POSITION FOUND: \t%s\t\t%s" % (each_vcf, directory))
                if not os.path.exists(directory):
                    try:
                        os.makedirs(directory)
                    except FileExistsError:
                        null = "null"
                shutil.copy(each_vcf, directory)
                ### ADD GROUP TO LIST
                group_calls.append(directory)
                
        # find mixed isolates if defining snp is ambigous
        for amb_position in list_amb:
            if amb_position in defining_snps:
                directory = defining_snps[amb_position]
                dict_amb.update({each_vcf + "\t" + directory:amb_position})
                ### ADD AMBIGIOUS CALL TO LIST
                group_calls.append("*" + directory + "-mix")
        # if -a or -e (non elites already deleted from the analysis) copy all vcfs to All_VCFs
        if args.all_vcf or args.elite:
            if not os.path.exists("All_VCFs"):
                os.makedirs("All_VCFs")
            shutil.move(each_vcf, "All_VCFs")
        else:
            try:
                os.remove(each_vcf)
            except FileNotFoundError:
                pass
        #print (dict_amb, group_calls, malformed)

    except ZeroDivisionError:
        os.remove(each_vcf)
        print ("ZeroDivisionError: corrupt VCF, removed %s " % each_vcf)
        mal = "ZeroDivisionError: corrupt VCF, removed %s " % each_vcf
        group_calls.append("error")
    except ValueError:
        os.remove(each_vcf)
        print ("ValueError: corrupt VCF, removed %s " % each_vcf)
        mal = "ValueError: corrupt VCF, removed %s " % each_vcf
        group_calls.append("error")
    except UnboundLocalError:
        os.remove(each_vcf)
        print ("UnboundLocalError: corrupt VCF, removed %s " % each_vcf)
        mal = "UnboundLocalError: corrupt VCF, removed %s " % each_vcf
        group_calls.append("error")
    except TypeError:
        os.remove(each_vcf)
        print ("TypeError: corrupt VCF, removed %s " % each_vcf)
        mal = "TypeError: corrupt VCF, removed %s " % each_vcf
        group_calls.append("error")
    except SyntaxError:
        os.remove(each_vcf)
        print ("SyntaxError: corrupt VCF, removed %s " % each_vcf)
        mal = "SyntaxError: corrupt VCF, removed %s " % each_vcf
        group_calls.append("error")
    except KeyError:
        os.remove(each_vcf)
        print ("KeyError: corrupt VCF, removed %s " % each_vcf)
        mal = "KeyError: corrupt VCF, removed %s " % each_vcf
        group_calls.append("error")
    except StopIteration:
        print ("StopIteration: %s" % each_vcf)
        mal = "KeyError: corrupt VCF, removed %s " % each_vcf
        group_calls.append("error")

    the_sample_name = group_calls[0:1]
    list_of_groups = sorted(group_calls[1:]) # order the groups
    for i in list_of_groups:
        the_sample_name.append(i) # a is group_calls
        group_calls = the_sample_name
    return dict_amb, group_calls, mal

# Group files, map pooled from script 2
def find_positions(filename):
    found_positions = {}
    vcf_reader = vcf.Reader(open(filename, 'r'))
    try:
        for record in vcf_reader:
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
    return found_positions

    def bruc_private_codes(upload_to):

    found = False
    if os.path.isfile("/Volumes/MB/Brucella/Brucella Logsheets/ALL_WGS.xlsx"):
        private_location = "/Volumes/MB/Brucella/Brucella Logsheets/ALL_WGS.xlsx"
        print("private_location:  %s" % private_location)
        found = True

    elif os.path.isfile("/fdrive/Brucella/Brucella Logsheets/ALL_WGS.xlsx"):
        private_location = "/fdrive/Brucella/Brucella Logsheets/ALL_WGS.xlsx"
        print("private_location:  %s" % private_location)
        found = True

    else:
        print("Path to Brucella genotyping codes not found")

    if found:
        wb_out = xlsxwriter.Workbook(upload_to + "/brucella/genotyping_codes.xlsx")
        ws_out = wb_out.add_worksheet()

        wb_in = xlrd.open_workbook(private_location)

        row = 0
        col = 0

        sheet_in = wb_in.sheet_by_index(1)
        for row_data in sheet_in.col(32):
            row_data = row_data.value
            row_data = re.sub("/", "_", row_data)
            row_data = re.sub("\.", "_", row_data)
            row_data = re.sub("\*", "_", row_data)
            row_data = re.sub("\?", "_", row_data)
            row_data = re.sub("\(", "_", row_data)
            row_data = re.sub("\)", "_", row_data)
            row_data = re.sub("\[", "_", row_data)
            row_data = re.sub("\]", "_", row_data)
            row_data = re.sub(" ", "_", row_data)
            row_data = re.sub("{", "_", row_data)
            row_data = re.sub("}", "_", row_data)
            row_data = re.sub("\'", "_", row_data)
            row_data = re.sub("-_", "_", row_data)
            row_data = re.sub("_-", "_", row_data)
            row_data = re.sub("--", "_", row_data)
            row_data = re.sub("_$", "", row_data)
            row_data = re.sub("-$", "", row_data)
            row_data = re.sub("\'", "", row_data)
            row_data = str(row_data)

            ws_out.write(row, col, row_data)
            row += 1

        wb_out.close()

table_location = outdir + directory + "-table.txt"
table=open(table_location, 'wt')

# write absolute positions to table
# order before adding to file to match with ordering of individual samples below
# all_positions is abs_pos:REF
all_positions=OrderedDict(sorted(all_positions.items()))

# Add the positions to the table
print ("reference_pos", end="\t", file=table)
for k, v in all_positions.items():
    print(k, end="\t", file=table)
print ("", file=table)

list_of_files = glob.glob('*vcf')

# for each vcf
all_map_qualities={}
for file_name in list_of_files:
    sample_map_qualities={}
    just_name = file_name.replace('.vcf', '')
    just_name = re.sub('\..*', '*', just_name) # if after the .vcf is removed there is stilll a "." in the name it is assumed the name did not get changed
    print(just_name, end="\t", file=table)
    # for each line in vcf
    vcf_reader = vcf.Reader(open(file_name, 'r'))
    sample_dict = {}
    for record in vcf_reader:
        record_position = str(record.CHROM) + "-" + str(record.POS)
        if record_position in all_positions:
            #print ("############, %s, %s" % (file_name, record_position))
            # NOT SURE THIS IS THE BEST PLACE TO CAPTURE MQ AVERAGE
            # MAY BE FASTER AFTER PARSIMONY SNPS ARE DECIDED, BUT THEN IT WILL REQUIRE OPENING THE FILES AGAIN.
            if str(record.ALT[0]) != "None" and str(record.INFO['MQ']) != "nan": #on rare occassions MQ gets called "NaN" thus passing a string when a number is expected when calculating average.
                #print ("getting map quality:    %s          %s      %s" % (record.INFO['MQ'], file_name, str(record.POS)))
                sample_map_qualities.update({record_position:record.INFO['MQ']})
            # ADD PARAMETERS HERE TO CHANGE WHAT'S EACH VCF REPRESENTS.
            # SNP IS REPRESENTED IN TABLE, NOW HOW WILL THE VCF REPRESENT THE CALLED POSITION
            # str(record.ALT[0]) != "None", which means a deletion as ALT
            # not record.FILTER, or rather PASSED.
            
            # check record.QUAL
            # In GATK VCFs "!= None" not used.
            if str(record.ALT[0]) != "None" and len(record.ALT[0]) == 1 and record.INFO['AC'][0] == 2 and record.QUAL > N_gatk_threshold:
                sample_dict.update({record_position:record.ALT[0]})
            elif str(record.ALT[0]) != "None" and len(record.ALT[0]) == 1 and record.INFO['AC'][0] == 1 and int(record.QUAL) > N_gatk_threshold:
                ref_alt = str(record.ALT[0]) + str(record.REF[0])
                if ref_alt == "AG":
                    sample_dict.update({record_position:"R"})
                elif ref_alt == "CT":
                    sample_dict.update({record_position:"Y"})
                elif ref_alt == "GC":
                    sample_dict.update({record_position:"S"})
                elif ref_alt == "AT":
                    sample_dict.update({record_position:"W"})
                elif ref_alt == "GT":
                    sample_dict.update({record_position:"K"})
                elif ref_alt == "AC":
                    sample_dict.update({record_position:"M"})
                elif ref_alt == "GA":
                    sample_dict.update({record_position:"R"})
                elif ref_alt == "TC":
                    sample_dict.update({record_position:"Y"})
                elif ref_alt == "CG":
                    sample_dict.update({record_position:"S"})
                elif ref_alt == "TA":
                    sample_dict.update({record_position:"W"})
                elif ref_alt == "TG":
                    sample_dict.update({record_position:"K"})
                elif ref_alt == "CA":
                    sample_dict.update({record_position:"M"})
                else:
                    sample_dict.update({record_position:"N"})
                # Poor calls
            elif str(record.ALT[0]) != "None" and int(record.QUAL) <= 50:
                sample_dict.update({record_position:record.REF[0]})
            elif str(record.ALT[0]) != "None" and int(record.QUAL) <= N_gatk_threshold:
                sample_dict.update({record_position:"N"})
            elif str(record.ALT[0]) != "None": #Insurance -- Will still report on a possible SNP even if missed with above statement
                sample_dict.update({record_position:str(record.REF[0])})
            elif str(record.ALT[0]) == "None":
                sample_dict.update({record_position:"-"})

    # After iterating through VCF combine dict to nested dict
    all_map_qualities.update({just_name: sample_map_qualities})

    # merge dictionaries and order
    merge_dict={}
    merge_dict.update(all_positions) #abs_pos:REF
    merge_dict.update(sample_dict) # abs_pos:ALT replacing all_positions, because keys must be unique
    merge_dict=OrderedDict(sorted(merge_dict.items())) #OrderedDict of ('abs_pos', ALT_else_REF), looks like a list of lists
    for k, v in merge_dict.items():
        #print ("k %s, v %s" % (k, v))
        print (str(v) + "\t", file=table, end="")
    print ("", file=table) # sample printed to file
table.close() #end of loop.  All files done

## Select parsimony informative SNPs
mytable = pd.read_csv(table_location, sep='\t')
# drop NaN rows and columns
mytable=mytable.dropna(axis=1)

# SELECT PARISOMONY INFORMATIVE SNPSs
# removes columns where all fields are the same
parsimony=mytable.loc[:, (mytable != mytable.iloc[0]).any()]
parsimony_positions=list(parsimony)
#write over table (table_location) containing all snps
parsimony.to_csv(table_location, sep="\t", index=False)

table=open(table_location, 'a')
# The reference calls are added after the parsimony positions are selected.
# added corresponding reference to parsimony table
print ("reference_call", end="\t", file=table)
#all_positions_list=list(all_positions)
try: #if there is only one file in the group exception is needed to return a value
    parsimony_positions.remove('reference_pos')
except ValueError:
    samples_in_fasta = []
    return(samples_in_fasta)

list_of_ref = []
for abs_pos in parsimony_positions:
    list_of_ref.append(all_positions.get(abs_pos))
string_of_ref = "\t".join(list_of_ref)
print(string_of_ref, file=table)
table.close()

samples_in_fasta = []
#Print out fasta alignment file from table
alignment_file= outdir + directory + ".fasta"
write_out=open(alignment_file, 'wt')
with open(table_location, 'rt') as f:
    count=0
    for line in f:
        if count > 0:
            line=re.sub('^', '>', line)
            line=line.replace('reference_call', 'root')
            line=line.replace('\t', '\n', 1)
            samples_in_fasta.append(line.split('\n')[0].replace('>', ''))
            line=line.replace('\t', '')
            print (line, end="", file=write_out)
        count = count + 1
write_out.close()

try: #if there are no SNP is the table
    mytable = pd.read_csv(table_location, sep='\t')
except:
    samples_in_fasta = []
    return(samples_in_fasta)

# move reference to top row
myref=mytable[-1:]
myother=mytable[:-1]
frames = [myref, myother]
mytable=pd.concat(frames)
mytable.to_csv(table_location, sep="\t", index=False)

print ("\n%s table dimensions: %s" % (directory, str(mytable.shape)))

print ("%s RAxML running..." % directory)
rooted_tree = outdir + directory + "-rooted.tre"
try:
    os.system("{} -s {} -n raxml -m GTRCATI -o root -p 12345 -T {} > /dev/null 2>&1" .format(sys_raxml, alignment_file, raxml_cpu))
except:
    write_out=open('RAXML_FAILED', 'w+')
    write_out.close()
    pass

def sort_table(table_location, ordered, out_org):
        mytable = pd.read_csv(table_location, sep='\t')
        #mytable=mytable.set_index('reference_pos')

        # order list is from tree file
        # gives order for samples to be listed in table to be phylogenetically correct
        ordered_list = []
        with open(ordered) as infile:
            for i in infile:
                i = i.rstrip()
                ordered_list.append(i)

        # Convert reference_pos-column to category and in set the ordered_list as categories hierarchy
        mytable.reference_pos = mytable.reference_pos.astype("category")
        mytable.reference_pos.cat.set_categories(ordered_list, inplace=True)
        mytable = mytable.sort_values(["reference_pos"]) # 'sort' changed to 'sort_values'

        # count number of SNPs in each column
        snp_per_column = []
        for column_header in mytable:
            count = 0
            column = mytable[column_header]
            # for each element in the column
            for element in column:
                if element != column[0]:
                    count = count + 1
            snp_per_column.append(count)
            #print ("the count is: %s" % count)
        row1 = pd.Series (snp_per_column, mytable.columns, name="snp_per_column")
        #row1 = row1.drop('reference_pos')

        # get the snp count per column
        # for each column in the table
        snp_from_top = []
        for column_header in mytable:
            count = 0
            column = mytable[column_header]
            # for each element in the column
            # skip the first element
            for element in column[1:]:
                if element == column[0]:
                    count = count + 1
                else:
                    break
            snp_from_top.append(count)
        row2 = pd.Series (snp_from_top, mytable.columns, name="snp_from_top")
        #row2 = row2.drop('reference_pos')

        mytable = mytable.append([row1])
        mytable = mytable.append([row2])
        
#In pandas=0.18.1 even this does not work:
#    abc = row1.to_frame()
#    abc = abc.T --> mytable.shape (5, 18), abc.shape (1, 18)
#    mytable.append(abc)
#Continue to get error: "*** ValueError: all the input arrays must have same number of dimensions"

        mytable = mytable.T
        mytable = mytable.sort_values(['snp_from_top', 'snp_per_column'], ascending=[True, False])
        mytable = mytable.T

        # remove snp_per_column and snp_from_top rows
        mytable = mytable[:-2]
        mytable.to_csv(out_org, sep='\t', index=False)

try:
    ordered_list_from_tree = outdir + directory + "-cleanedAlignment.txt"
    write_out=open(ordered_list_from_tree, 'w+')
    print ("reference_pos", file=write_out)
    print ("reference_call", file=write_out)
    if os.path.isfile("RAxML_bestTree.raxml"):
        with open("RAxML_bestTree.raxml", 'rt') as f:
            for line in f:
                line=re.sub('[:,]', '\n', line)
                line=re.sub('[)(]', '', line)
                line=re.sub('[0-9].*\.[0-9].*\n', '', line)
                line=re.sub('root\n', '', line)
                write_out.write(line)
        best_raxml_tre = directory + "-RAxML-bestTree.tre"
        os.rename("RAxML_bestTree.raxml", best_raxml_tre)
        write_out.close()

    best_raxml_svg = directory + "-RAxML-bestTree.svg"
    best_raxml_pdf = directory + "-RAxML-bestTree.pdf"
    
    try:
        os.system("cat {} | nw_display -s -S -w 1300 -t -v 30 -i 'opacity:0' -b 'opacity:0' -l 'font-size:14;font-family:serif;font-style:italic' -d 'stroke-width:1;stroke:blue' - > {}" .format(best_raxml_tre, best_raxml_svg)) #-s produces svg, -S suppress scale bar, -w to set the number of columns available for display, -t tab format, -v vertical spacing, -i inner node label, -b branch style
        svg2pdf(url=best_raxml_svg, write_to=best_raxml_pdf)
    except:
        pass
    
    out_org = outdir + directory + "-organized-table.txt"

    sort_table(table_location, ordered_list_from_tree, out_org) #function

    print ("%s Getting map quality..." % directory)
    average=lambda x: x.mean()
    all_map_qualities=pd.DataFrame(all_map_qualities)
    #ave_mq = Type: Series
    ave_mq = all_map_qualities.apply(average, axis=1)
    ave_mq = ave_mq.astype(int)
    ave_mq.to_csv('outfile.txt', sep='\t') # write to csv

    write_out=open('map_quality.txt', 'w+')
    print ('reference_pos\tmap-quality', file=write_out)
    with open('outfile.txt', 'rt') as f:
        for line in f:
            write_out.write(line)
    write_out.close()
    
    #seemed pooling did not like a function with no parameters given
    quality = pd.read_csv('map_quality.txt', sep='\t')

    mytable = pd.read_csv(table_location, sep='\t')
    mytable=mytable.set_index('reference_pos')

    # order list is from tree file
    # gives order for samples to be listed in table to be phylogenetically correct
    ordered_list = []
    with open(ordered_list_from_tree) as infile:
        for i in infile:
            i = i.rstrip()
            ordered_list.append(i)
    # sinces this is set as the mytable index do not include in ordering
    ordered_list.remove('reference_pos')

    # reorder table based on order of list
    mytable = mytable.reindex(ordered_list)
    mytable.to_csv(table_location, sep='\t')

    out_sort=str(os.getcwd()) + "/" + directory + "-sorted-table.txt" #sorted
    mytable_sort = pd.read_csv(table_location, sep='\t') #sorted
    mytable_sort = mytable_sort.set_index('reference_pos') #sorted
    mytable_sort = mytable_sort.transpose() #sort
    mytable_sort.to_csv(out_sort, sep='\t', index_label='reference_pos') #sort

    out_org=str(os.getcwd()) + "/" + directory + "-organized-table.txt" #org
    mytable = pd.read_csv(out_org, sep='\t') #org
    mytable = mytable.set_index('reference_pos') #org
    mytable = mytable.transpose() #org
    mytable.to_csv(out_org, sep='\t', index_label='reference_pos') #org

    if mygbk and not args.no_annotation:
        dict_annotation = get_annotations_table(parsimony_positions)
        write_out=open('annotations.txt', 'w+')
        print ('reference_pos\tannotations', file=write_out)
        for k, v in dict_annotation.items():
            print ('%s\t%s' % (k, v), file=write_out)
        write_out.close()
    
        print ("%s gbk is present, getting annotation...\n" % directory)
        annotations = pd.read_csv('annotations.txt', sep='\t') #sort
        mytable_sort = pd.read_csv(out_sort, sep='\t') #sort
        mytable_sort = mytable_sort.merge(quality, on='reference_pos', how='inner') #sort
        mytable_sort = mytable_sort.merge(annotations, on='reference_pos', how='inner') #sort
        mytable_sort = mytable_sort.set_index('reference_pos') #sort
        mytable_sort = mytable_sort.transpose() #sort
        mytable_sort.to_csv(out_sort, sep='\t', index_label='reference_pos') #sort

        #annotations = pd.read_csv('annotations.txt', sep='\t') #org
        mytable = pd.read_csv(out_org, sep='\t') #org
        mytable = mytable.merge(quality, on='reference_pos', how='inner') #org
        mytable = mytable.merge(annotations, on='reference_pos', how='inner') #org
        mytable = mytable.set_index('reference_pos') #org
        mytable = mytable.transpose() #org
        mytable.to_csv(out_org, sep='\t', index_label='reference_pos') #org

    else:
        print ("No gbk file or no table to annotate")
        mytable_sort = pd.read_csv(out_sort, sep='\t') #sort
        mytable_sort = mytable_sort.merge(quality, on='reference_pos', how='inner') #sort
        mytable_sort = mytable_sort.set_index('reference_pos') #sort
        mytable_sort = mytable_sort.transpose() #sort
        mytable_sort.to_csv(out_sort, sep='\t', index_label='reference_pos') #sort
        # add when no annotation
        with open(out_sort, 'rt') as f:
            line=f.readline()
        f.close()
        column_count = line.count('\t') #sort
        column_count = column_count - 1 #sort
        #print ("column_count: %s" % column_count)
        with open(out_sort, 'at') as f:
            print ("no_annotation", end = '', file=f)
            print ('\t' * column_count, file=f)
        f.close()

        print ("No gbk file or no table to annotate")
        mytable = pd.read_csv(out_org, sep='\t') #org
        mytable = mytable.merge(quality, on='reference_pos', how='inner') #org
        mytable = mytable.set_index('reference_pos') #org
        mytable = mytable.transpose() #org
        mytable.to_csv(out_org, sep='\t', index_label='reference_pos') #org
        # add when no annotation
        with open(out_org, 'rt') as f:
            line=f.readline()
        f.close()
        column_count = line.count('\t')
        column_count = column_count - 1
        #print ("column_count: %s" % column_count)
        with open(out_org, 'at') as f:
            print ("no_annotation", end = '', file=f)
            print ('\t' * column_count, file=f)
        f.close()

    def excelwriter(filename):
        orginal_name=filename
        filename = filename.replace(".txt",".xlsx")
        wb = xlsxwriter.Workbook(filename)
        ws = wb.add_worksheet("Sheet1")
        with open(orginal_name,'r') as csvfile:
            table = csv.reader(csvfile, delimiter='\t')
            i = 0
            for row in table:
                ws.write_row(i, 0, row)
                i += 1

        col = len(row)
        col = col + 1
        #print (i, "x", col)

        formatA = wb.add_format({'bg_color':'#58FA82'})
        formatG = wb.add_format({'bg_color':'#F7FE2E'})
        formatC = wb.add_format({'bg_color':'#0000FF'})
        formatT = wb.add_format({'bg_color':'#FF0000'})
        formatnormal = wb.add_format({'bg_color':'#FDFEFE'})
        formatlowqual = wb.add_format({'font_color':'#C70039', 'bg_color':'#E2CFDD'})
        formathighqual = wb.add_format({'font_color':'#000000', 'bg_color':'#FDFEFE'})
        formatambigous = wb.add_format({'font_color':'#C70039', 'bg_color':'#E2CFDD'})
        formatN = wb.add_format({'bg_color':'#E2CFDD'})

        ws.conditional_format(i-2,1,i-2,col-2, {'type':'text',
                                'criteria':'containing',
                                'value':60,
                                'format':formathighqual})
        ws.conditional_format(i-2,1,i-2,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':59,
                            'format':formathighqual})
        ws.conditional_format(i-2,1,i-2,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':58,
                            'format':formathighqual})
        ws.conditional_format(i-2,1,i-2,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':57,
                            'format':formathighqual})
        ws.conditional_format(i-2,1,i-2,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':56,
                            'format':formathighqual})
        ws.conditional_format(i-2,1,i-2,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':55,
                            'format':formathighqual})
        ws.conditional_format(i-2,1,i-2,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':54,
                            'format':formathighqual})
        ws.conditional_format(i-2,1,i-2,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':53,
                            'format':formathighqual})
        ws.conditional_format(i-2,1,i-2,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':52,
                            'format':formathighqual})
        ws.conditional_format(i-2,1,i-2,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':51,
                            'format':formathighqual})
        ws.conditional_format(i-2,1,i-2,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':50,
                            'format':formathighqual})
        ws.conditional_format(i-2,1,i-2,col-2, {'type':'text',
                            'criteria':'not containing',
                            'value':100,
                            'format':formatlowqual})

        ws.conditional_format(2,1,i-3,col-2, {'type':'cell',
                            'criteria':'==',
                            'value':'B$2',
                            'format':formatnormal})
        ws.conditional_format(2,1,i-3,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':'A',
                            'format':formatA})
        ws.conditional_format(2,1,i-3,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':'G',
                            'format':formatG})
        ws.conditional_format(2,1,i-3,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':'C',
                            'format':formatC})
        ws.conditional_format(2,1,i-3,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':'T',
                            'format':formatT})
        ws.conditional_format(2,1,i-3,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':'S',
                            'format':formatambigous})
        ws.conditional_format(2,1,i-3,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':'Y',
                            'format':formatambigous})
        ws.conditional_format(2,1,i-3,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':'R',
                            'format':formatambigous})
        ws.conditional_format(2,1,i-3,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':'W',
                            'format':formatambigous})
        ws.conditional_format(2,1,i-3,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':'K',
                            'format':formatambigous})
        ws.conditional_format(2,1,i-3,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':'M',
                            'format':formatambigous})
        ws.conditional_format(2,1,i-3,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':'N',
                            'format':formatN})
        ws.conditional_format(2,1,i-3,col-2, {'type':'text',
                            'criteria':'containing',
                            'value':'-',
                            'format':formatN})

        ws.set_column(0, 0, 30)
        ws.set_column(1, col-2, 2)
        ws.freeze_panes(2, 1)
        format_rotation = wb.add_format({'rotation':'90'})
        ws.set_row(0, 140, format_rotation)
        formatannotation = wb.add_format({'font_color':'#0A028C', 'rotation':'-90', 'align':'top'})
        #set last row
        ws.set_row(i-1, 400, formatannotation)

        wb.close()

    excelwriter(out_sort) #***FUNCTION CALL #sort
    excelwriter(out_org) #***FUNCTION CALL #org

    for r in glob.glob('*vcf'):
        os.remove(r)

except ValueError:
    print ("##### ValueError: %s #####" % file_name)
    return

try:
    os.remove(ordered_list_from_tree)
    os.remove('map_quality.txt')
    if mygbk:
        os.remove("annotations.txt")
    os.remove("outfile.txt")
    os.remove(out_sort)
    os.remove(out_org) # organized.txt table
    os.remove(table_location) # unorganized table
    os.remove('RAxML_info.raxml')
    os.remove('RAxML_log.raxml')
    os.remove('RAxML_parsimonyTree.raxml')
    os.remove('RAxML_result.raxml')
    os.remove(directory + '.fasta.reduced')

except FileNotFoundError:
    pass

### PANDA NOTES ###
# get the index: mytable.index
# get columns: mytable.columns
# get a column: mytable.AF2122_NC002945_105651, shows index (sample names)
# get a row: mytable.ix['reference'], shows columns (positions and SNPs)
# values: mytable.values, SNPs - series
# strip off the bottom row: mytable[:-1]
# get the bottom row: mytable[-1:]

with open(directory + "samples_in_fasta.json", 'w') as outfile:
    json.dump(samples_in_fasta, outfile)

return(samples_in_fasta)
