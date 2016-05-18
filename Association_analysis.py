#!/user/bin/env python
# Filename = Association_analysis.py
# -*- coding=Unicode -*- 

def List_disease(folderPath):     # The function can only process HMDB metabolites database. The output file is sorted by diseases.
    import re, bs4, os

    if folderPath[-1] != '/':
        folderPath += '/'
    filename_list = os.listdir(folderPath)
    data = []
    diseases_names = []
    for files in filename_list:
        print "loading%s..." % files
        f = open(folderPath + files, 'r')
        soup = bs4.BeautifulSoup(f, "html.parser")
        if not soup.abnormal_concentrations.string:
            meta_name = re.findall(r'<name>([^<]+)</name>', str(soup.metabolite))[0]
            if soup.tissue_locations.string:
                Tissue = 'Not available'
            else:
                Tissue = ', '.join(re.findall(r'<tissue>([^<]+)</tissue>', str(soup.tissue_locations)))
            if soup.pathways.string:
                Pathway = 'Not available'
            else:
                namesid = list(re.findall(r'<name>([^<]+)</name>\n<smpdb_id>([0-9a-zA-Z]*)</smpdb_id>\n<kegg_map_id>([0-9a-zA-Z]*)</kegg_map_id>', str(soup.pathways)))
                nameids = []
                for item in namesid:
                    nameids.append('%s(%s,%s)'% (item[0],item[1],item[2]))
                Pathway = ', '.join(nameids)

            Nor_para = re.findall(r'<biofluid>([^<]+)</biofluid>\n<concentration_value>[ &;a-zA-Z]*([.0-9]+)[ ().+-/0-9]*</concentration_value>\n<concentration_units>[^<]+</concentration_units>\n<subject_age>([^<]+)</subject_age>\n<subject_sex>([^<]+)</subject_sex>\n<subject_condition>([^<]+)</subject_condition>', str(soup.normal_concentrations))
            Nor_para = map(list, Nor_para)      # Turn tuple into list.
            Nor_para = sorted(Nor_para, key=lambda Nor_para : (Nor_para[4], Nor_para[3], Nor_para[2], Nor_para[0]))
            Abn_para = re.findall(r'<biofluid>([^<]+)</biofluid>\n<concentration_value>[ &;a-zA-Z]*([.0-9]+)[ .+-/0-9]*</concentration_value>\n<concentration_units>[^<]+</concentration_units>\n<patient_age>([^<]+)</patient_age>\n<patient_sex>([^<]+)</patient_sex>\n<patient_information>([^<]+)</patient_information>', str(soup.abnormal_concentrations))
            Abn_para = map(list, Abn_para)
            Abn_para = sorted(Abn_para, key=lambda Abn_para : (Abn_para[4], Abn_para[3], Abn_para[2], Abn_para[0]))

            tempData = {}
            for i in xrange(0, len(Abn_para)):
                if Abn_para[i][4].strip() not in tempData:
                    tempData[Abn_para[i][4].strip()] = [0,0]
                for No in Nor_para:
                    if No[0] == Abn_para[i][0] and No[2] == Abn_para[i][2] and No[3] == Abn_para[i][3]:
                        if float(Abn_para[i][1]) > float(No[1]):
                            tempData[Abn_para[i][4].strip()][0] += 1
                        if float(Abn_para[i][1]) < float(No[1]):
                            tempData[Abn_para[i][4].strip()][1] += 1
            for item in tempData:
                if tempData[item][0] > 0:
                    if tempData[item][1] == 0:
                        tempData[item] = '+'
                    else:
                        tempData[item] = '%s' %  (round((float(tempData[item][0]-tempData[item][1])/float(tempData[item][0]+tempData[item][1])),4))
                elif tempData[item][1] > 0:
                    tempData[item] = '-'
                data.append([item, meta_name.strip(), Tissue.strip(), Pathway.strip(), tempData[item]])
                if item not in diseases_names:
                    diseases_names.append(item)
        f.close()

    diseases_names = sorted(diseases_names, key=lambda x : x)
    
    data = sorted(data, key=lambda x : x[0])
    for i in xrange(len(diseases_names)):
        print i, diseases_names[i]

    judge = False
    while judge == False:
        userInput = raw_input('Input the number index of what diseases you need from those diseases names above\nSeperate diseases with comma if they belong to the same category, and seperate categories with "/"\nIf input "all" you can get all the data\n')
        if not re.search(r'[^0-9/,]',userInput) or userInput == 'all':
            judge = True
    paras = []
    if not os.path.exists('step1'):
        os.mkdir('step1')
    if userInput.lower() == 'all':
        f = open('step1/Alldiseases.dat', 'w')
        for content in data:
            f.write('%s\t%s\t%s\t%s\t%s\n' % (content[0], content[1], content[2], content[3], content[4]))
    else:
        for item in userInput.strip().split('/'):
            if item != '':
                tempList = item.split(',')
                for thing in tempList:
                    thing = thing.strip()
                paras.append(tempList)

        for i in xrange(len(paras)):
            exec("f = open('step1/Dis%s.dat', 'w')" % str(i+1))      # i+1 is for convenience.
            string = '#Disease\tMetabolite\tTissue location\tPathway\tChange of Concentration\\n'
            exec("f.write('%s')" % string)

            name = []
            for index in paras[i]:
                name.append(diseases_names[int(index)])
            for content in data:
                if content[0] in name:
                    f.write('%s\t%s\t%s\t%s\t%s\n' % (content[0], content[1], content[2], content[3], content[4]))
            f.close()

def transform(listName, mainparaloc=3, minorparasloc=(0,1,4), keyword='map'):   # This function supplements List_pathway so don't change it or use it solely.
    import re
    Dict = {}
    for item in listName:
        paras = item.strip().split('\t')
        mainpara = paras[mainparaloc]
        if keyword not in mainpara and 'Not available' in mainpara:
            mainkeys = ['Not available']
        elif keyword in mainpara:
            mainkeys = re.findall('%s(\d+)' % keyword, mainpara)
        else:
            continue

        for key in mainkeys:
            if key not in Dict:
                Dict[key] = '{%s,%s,%s}'% eval('(paras[%d],paras[%d],paras[%d])' % minorparasloc)
            else:
                Dict[key]+='{%s,%s,%s}'% eval('(paras[%d],paras[%d],paras[%d])' % minorparasloc)
    return Dict

def List_pathway(folderPath, Not_available = True, All_Tissues = True): 
    import os,re
    if not os.path.exists('step2'):
        os.mkdir('step2')
    if folderPath[-1] != '/':
        folderPath += '/'
    filename_list = os.listdir(folderPath)
    locList = []
    for files in filename_list:
        f = open(folderPath + files, 'r')
        for line in f:
            if not re.match('#', line):
                for item in line.split('\t')[2].strip().split(', '):
                    if item not in locList:
                        locList.append(item)
        f.close()
    if 'All Tissues' in locList:
        locList.remove('All Tissues')
    if 'Not available' in locList:
        locList.remove('Not available')
    locList = sorted(locList, key=lambda loc: loc[0])
    for i in xrange(len(locList)):
        print i, locList[i]

    Dict = {}
    DictN = {}
    userInput = raw_input('Input the index of what tissue_locations you need from those location names above\nSeperate locations with comma if they belong to the same category, and seperate categories with "/"\n')
    locs = userInput.strip().split('/')
    for i in xrange(len(locs)):
        Dict[i] = []
        for locindex in locs[i].split(','):
            Dict[i].append(locList[int(locindex.strip())])
        if Not_available == True:
            Dict[i].append('Not available')
        if All_Tissues == True:
            Dict[i].append('All Tissues')

    for files in filename_list:
        for i in xrange(len(locs)):
            if not os.path.exists('step2/%s' %str(i)):
                os.mkdir('step2/%s' %str(i))
            exec("f%s = open('step2/%s/%s_Tis%s', 'w')" % (i, i, files.split('/')[-1].split('.')[0], i))
            exec("L%s = []" % i)
            DictN[i] = 'L%s.append(line)' % i

        f = open(folderPath + files, 'r')
        for line in f:
            if not re.match('#', line):
                for item in Dict:
                    for i in Dict[item]:
                        if i in line:
                            exec(DictN[item])
                            break

        for i in xrange(len(locs)):
            temp = sorted(eval('transform(L%s)' % i).iteritems(), key=lambda d:d[0], reverse = True)
            for item in temp:
                string = '%s\t%s\n' % (item[0], item[1])
                exec('f%s.write(string)' % i)
            exec("f%s.close()" % i)
        f.close()

###############################################################################################################

def integrate_updown(folderPath):               # It only works for symbols file of NCBI GEO dataset.
    import re, os,time
    if not os.path.exists('step3'):
        os.mkdir('step3')
    if folderPath[-1] != '/':
        folderPath += '/'
    filename_list = os.listdir(folderPath)
    res = {}
    foldcol = 0
    f = open('step3/Res_%s' % folderPath.split('/')[-2], 'w')
    for files in filename_list:
        F = open(folderPath + files, 'r')
        paras = F.readline().strip().split('\t')
        for i in xrange(len(paras)):
            if paras[i] == 'fold':
                foldcol = i

        for line in F:
            if not re.match('Proset', line):
                para = line.strip().split('\t')
                gene = para[1].split(':')[-1]
                fold = para[foldcol]

                if gene not in res:
                    if '-' in fold:
                        res[gene] = -1
                    else:
                        res[gene] = 1
                else:
                    if '-' in fold:
                        res[gene] -= 1 
                    else:
                        res[gene] += 1
        F.close()

    for item in res:
        if res[item] > 0:
            f.write('%s\t+\n' % item)
        elif res[item] < 0:
            f.write('%s\t-\n' % item)
        else:
            f.write('%s\t/\n' % item)
    f.close()

def match(metaFile, tranFile, updownFile):   # metaFile and updownFile are both outputs of the functions above, while tranFile need to be produced with DAVID tools.
    import re,os
    metaExist = tranExist = False
    if metaFile != '':
        F1 = open(metaFile, 'r')
        metaExist = True
    if tranFile != '' and updownFile != '':
        F2 = open(tranFile, 'r')
        F3 = open(updownFile, 'r')
        tranExist = True
    if not os.path.exists('step4'):
        os.mkdir('step4')
    if metaExist == True:
        if tranExist == True:
            f = open('step4/metran_%s.dat' % metaFile.split('/')[-1].split('.')[0], 'w')
        else:
            f = open('step4/meta_%s.dat' % metaFile.split('/')[-1].split('.')[0], 'w')
    elif tranExist == True:
        f = open('step4/tran_%s.dat' % tranFile.split('/')[-1].split('.')[0], 'w')
    f.write('#pathway\tmetabolome\tgene\n')

    if tranExist == True:
        F3Dict = {}
        for F3line in F3:
            if not re.match('#', F3line):
                paras = F3line.split()
                F3Dict[paras[0].upper()] = paras[1]
        F2Dict = {}
        for F2line in F2:
            if not re.match('Category', F2line):
                paras = F2line.split('\t')
                tranKegg = re.findall('(\d+):', paras[1])[0]
                genes = paras[5].split(', ')
                for i in xrange(len(genes)):
                    if genes[i].upper() in F3Dict:
                        genes[i] = '{%s,%s}'%(genes[i].upper(),F3Dict[genes[i]])
                    else:
                        genes[i] = '{%s,/}' % genes[i].upper()
                F2Dict[tranKegg] = ''.join(genes[i] for i in xrange(len(genes)))

        if metaExist == True:
            F1Dict = {}
            for F1line in F1:
                if (not re.match('#', F1line)) and (not re.match('Not available', F1line)):
                    metaKegg = F1line.split('\t')[0]
                    F1Dict[metaKegg] = F1line.strip()
                elif re.match('Not available', F1line):
                    f.write(F1line)
            for item in F2Dict:
                if item in F1Dict:
                    f.write('%s\t%s\n' % (F1Dict[item], F2Dict[item]))
                else:
                    f.write('%s\t/\t%s\n' % (item, F2Dict[item]))
            F1.close()
        else:
            for item in F2Dict:
                f.write('%s\t/\t%s\n' % (item, F2Dict[item]))
        F2.close()
        F3.close()
    else:
        if metaExist == True:
            for F1line in F1:
                metaKegg = F1line.split('\t')[0]                
                if metaKegg == 'Not available':
                    f.write(F1line)
                else:
                    f.write('%s\t/\n' % F1line.strip())
            F1.close()
    f.close()


def matchs(metaFolder, tranFile, updownFile):
    import os
    if metaFolder == '':
        match('', tranFile, updownFile)
    else:
        if metaFolder[-1] != '/':
            metaFolder += '/'
        filename_list = os.listdir(metaFolder)
        for files in filename_list:
            metaFile = metaFolder + files
            match(metaFile, tranFile, updownFile)

# Here we get the final file with metabolites, genes and their relationship with the specified diseases in. 
###########################################################################################################

def get_metacode(metaName): # Supplemental function to LdFile_DlPic.
    import urllib, requests, bs4, re
    url = 'http://www.kegg.jp/dbget-bin/www_bfind_sub?mode=bfind&mode=bfind&max_hit=1000&locale=en&serv=kegg&dbkey=compound&keywords=%s&page=1' % urllib.quote_plus(metaName)
    res = requests.get(url) #download the web page
    res.raise_for_status()  #get the connection status
    soup = bs4.BeautifulSoup(res.text, "html.parser")
    divs = soup.find_all('div')
    result = []
    metaName = metaName.replace('(', '\(').replace(')', '\)')
    try:
        for item in divs:
            if len(item) == 6:
                if re.search(r' %s[ ;<]' % metaName, str(item.contents)):
                    result.append(re.findall(r'cpd:(C[0-9]+)', str(item.contents))[0])
        return result
    except:
        return result

def get_mapandstr(linecontent):   # Supplemental function to LdFile_DlPic. map means the num of map and the str means the content to give.
    import re
    mapandstr = []
    paras = linecontent.split('\t')
    mapandstr.append(paras[0])
    string = ''
    
    if len(paras) == 3:
        metas = re.findall(r"{(.*?)}", paras[1])
        if len(metas) > 0:
            for item in metas:
                if item.split(',')[-1] != ' 0]':
                    commaPos=list(((pos) for pos,val in enumerate(item) if val == ','))
                    meta = [item[:commaPos[0]],item[commaPos[0]+1:commaPos[-1]],item[commaPos[-1]+1:]]
                    metaName = meta[1]
                    metaRel = meta[2]
                    if metaRel == '-':
                        metaColor = 'cyan'
                    elif metaRel == '+':
                        metaColor = 'red'
                    else:
                        if float(metaRel)<0:
                            metaColor = 'cyan'
                        elif float(metaRel)>0:
                            metaColor = 'red'
                        else:
                            metaColor = 'yellow'
                    for code in get_metacode(metaName):
                        string += '%s %s\n' % (code, metaColor)
        trans = re.findall(r"{(.*?)}", paras[2])
        if len(trans) > 0:
            for item in trans:
                tran = item.split(',')
                tranName = tran[0]
                tranRel = tran[1]
                if tranRel == '-':
                    tranColor = 'cyan'
                elif tranRel == '+':
                    tranColor = 'red'
                else:
                    tranColor = 'yellow'
                string += '%s %s\n' % (tranName, tranColor)
        mapandstr.append(string)
        return mapandstr

    elif len(paras) == 2:
        metas = re.findall(r"{(.*?)}", paras[1])
        for item in metas:
            if item.split(',')[-1] != ' 0]':
                commaPos=list(((pos) for pos,val in enumerate(item) if val == ','))
                meta = [item[:commaPos[0]],item[commaPos[0]+1:commaPos[-1]],item[commaPos[-1]+1:]]
                metaName = meta[1]
                metaRel = meta[2]
                if metaRel == '-':
                    metaColor = 'cyan'
                elif metaRel == '+':
                    metaColor = 'red'
                else:
                    if float(metaRel)<0:
                        metaColor = 'cyan'
                    elif float(metaRel)>0:
                        metaColor = 'red'
                    else:
                        metaColor = 'yellow'
                for code in get_metacode(metaName):
                    string += '%s %s\n' % (code, metaColor)
        mapandstr.append(string)
        return mapandstr

def LdFile_DlMap(filePath, org_name='hsa'):   # Load the final file and download the pathway map.
    import re,time,urllib,os,sys
    from selenium import webdriver
    from selenium.webdriver.common.keys import Keys

    current_path = os.path.split(os.path.realpath(__file__))[0]
    sys.path.append(current_path)

    F = open(filePath, 'r')
    NAExist = False
    for line in F:
        if re.match('Not available', line):
            NAcontent = get_mapandstr(line)[1]    # First, load the content of Not available metas.
            NAExist = True
            F.seek(0)
            break
    F.seek(0)
    for line in F:
        if (not re.match('#', line)) and (not re.match('Not available', line)):
            (mapkey, content) = get_mapandstr(line)
            if NAExist == True:
                content += NAcontent
            content = content.replace('  red\n','').replace('  cyan\n','').replace('  yellow\n','')    # Delete those useless columns.
            url = 'http://www.genome.jp/kegg-bin/show_pathway?org_name=%s&mapno=%s&mapscale=1.0&show_description=hide' % (org_name,mapkey)
            driver = webdriver.Chrome()
            driver.get(url)
            now_handle = driver.current_window_handle

            html = urllib.urlopen(url)
            scode = html.read()
            if re.findall('does not exist</h3>',scode) != []:
                print '------------\nThe %s doesn\'t exist\n----------------' % org_name+mapkey
                continue

            driver.find_element_by_link_text("User data mapping").click()
            handles = driver.window_handles
            for handle in handles:
                if handle != now_handle:
                    driver.switch_to_window(handle)
                    elem = driver.find_element_by_css_selector('body > table > tbody > tr:nth-child(2) > td > table > tbody > tr > td:nth-child(1) > textarea')
                    elem.send_keys(content)
                    driver.find_element_by_css_selector('body > table > tbody > tr:nth-child(2) > td > input[type="submit"]:nth-child(4)').click()
            driver.switch_to_window(now_handle)
            img = driver.find_element_by_css_selector('body > img')
            src = img.get_attribute('src')
            urllib.urlretrieve(src, "%smap%s.png" % (filePath.split('/')[-1].split('.')[0],mapkey))
            driver.quit()
            content = ''

def LdFile_DlMaps(folderPath, org_name):
    import os
    if folderPath[-1] != '/':
        folderPath += '/'
    filename_list = os.listdir(folderPath)
    for files in filename_list:
        filePath = folderPath + files
        LdFile_DlMap(filePath, org_name)

def main():
    import os
    print '----------------------------------------------------------------------------------\n'
    print '1:List_disease, the function integrates HMDB dataset and lists all metabolites by their linked diseases\n'
    print '2:List_pathway, the function relists metabolites files by the pathways\n'
    print '3:integrate_updown, the function integrate the proportion between gene expression and diseases\n'
    print '4:match, the function integrate the metabolome and transcriptome files\n'
    print '5:LdFile_DlMap, the function load the result file and download the KEGG maps\n'
    print '6:Exit'
    func = raw_input('Input the number index of the function you need\n')

    if func == '1':
        folderPath0 = raw_input('Please input the path of folder of HMDB metabolites files which are extracted.\n')
        folderPath0 = unicode(folderPath0,"gbk").replace('\\', '/')
        List_disease(folderPath0)

    elif func == '2':
        folderPath1 = raw_input('You may directly press "Enter" to use the step1 folderPath, or input the path of another folder\n')
        folderPath1 = unicode(folderPath1,"gbk").replace('\\', '/')
        if folderPath1 == '':
            folderPath1 = 'step1'
        print 'Do you want those metabolites, the tissue locations of which are not available, to be counted in all pathway?\n'
        choice = raw_input('Input Y or N, and the former is default\n')
        if choice.lower() == 'y':
            Not_available = True
        elif choice.lower() == 'n':
            Not_available = False
        List_pathway(folderPath1, Not_available, All_Tissues = True)
        
    elif func == '3':
        folderPath2 = raw_input('Please input the path of folder of NCBI GEO proportion files.\n')
        folderPath2 = unicode(folderPath2,"gbk").replace('\\', '/')
        integrate_updown(folderPath2)
        
    elif func == '4':
        folderPath3 = raw_input('Please input the path of folder of metabolites result from the same location.\n')
        folderPath3 = unicode(folderPath3,"gbk").replace('\\', '/')
        tranFile = raw_input('Please input the path of GO/KEGG analysis result file.\n')
        tranFile = unicode(tranFile,"gbk").replace('\\', '/')
        updownFile = raw_input('Please input the path of corresponding function3 result file.\n')
        updownFile = unicode(updownFile,"gbk").replace('\\', '/')
        matchs(folderPath3, tranFile, updownFile)
        
    elif func == '5':
        folderPath4 = raw_input('You may directly press "Enter" to use the step4 folderPath, or input the path of another folder\n')
        if folderPath4 == '':
            folderPath4 = 'step4'
        else:
            folderPath4 = unicode(folderPath4,"gbk").replace('\\', '/')
        org_name = raw_input('You may directly press "Enter" to set Human as organism, or input the codes after the "?org_name=" and before the "&map" of the website\n')
        if org_name == '':
            org_name = 'hsa'
        LdFile_DlMaps(folderPath4, org_name)
        
    elif func == '6':
        exit()

while True:
    main()
