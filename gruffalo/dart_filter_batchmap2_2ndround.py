#!/usr/bin/env python3
    
#submit from /projects/hop/dartseq
#ProgDir=/home/hajduk/git_repos/scripts/hop_genomics/dart_filtering
#sbatch $ProgDir/dart_filter_batchmap2_2ndround.py   

#### final version used for linkage mapping 24.11.2021

import csv

def check_all_good(row,startindex):
    si=startindex
    if(row[si]==row[si+1] and row[si+1]==row[si+2]):
        return row[si]
    else:
        return "error"

def check_3_good(row,startindex):
    si=startindex
    counter=[0,0,0]
    counterconf=0
    for i in [si,si+1,si+2,si+3]:
        if(row[i]=="-"):
            counterconf+=1
        else:
            counter[int(row[i])]+=1
            if(counter[int(row[i])]==3): 
                return row[i]
    if(counterconf>=3):  
        return "--"    
    if(counterconf==2 and row[i]==2):
	    return "-"    
    return "error"        

def get_segtype(typeA,typeB):

#standard seg types  
    if(typeA=="2" and typeB=="2"):
        return 'B3.7'
    if(typeA=="2" and typeB=="0"):
        return 'D1.10'
    if(typeA=="2" and typeB=="1"):
        return 'D1.10b'
    if(typeA=="0" and typeB=="2"):
        return 'D2.15'
    if(typeA=="1" and typeB=="2"):
        return 'D2.15b'

#het x two null
    if(typeA=="2" and typeB=="--"):
        return 'D1.11'
    if(typeA=="--" and typeB=="2"):
        return 'D2.16'

# het x one null allele
    if(typeA=="2" and typeB=="-"):
        return 'B1.5'
    if(typeA=="-" and typeB=="2"):
        return 'B2.6'

# homo x one null allele
    if(typeA=="0" and typeB=="-"):
        return 'D2.17'
    if(typeA=="1" and typeB=="-"):
        return 'D2.17b'
    if(typeA=="-" and typeB=="0"):
        return 'D1.12'
    if(typeA=="-" and typeB=="1"):
        return 'D1.12b'

#accepted missing genotypes    
    if(typeA=="0" and typeB=="1"):
        return 'A4'
    if(typeA=="1" and typeB=="0"):
        return 'A4b'
    if(typeA=="-" and typeB=="--"):
        return 'D1.13'
    if(typeA=="--" and typeB=="-"):
        return 'D2.18'
    if(typeA=="0" or typeA=="1") and (typeB=="--"):
        return 'D1.13x'
    if(typeA=="--") and (typeB=="0" or typeB=="1"):
        return 'D2.18x'
    return "undef"              #after this argument is passed it just returns everything as undef

def convert_segformat_tooutput(cell,rowtype,output):

#standard segtypes 
    if(rowtype=="D1.10"):
        values=['a,','-,','ab,']
        if(cell=="-"):
            output.write("-,")
        else:
            output.write(values[int(cell)])
    if(rowtype=="D1.10b"):
        values=['-,','a,','ab,']  
        if(cell=="-"):
            output.write("-,")
        else:
            output.write(values[int(cell)])
    if(rowtype=="D2.15"):
        values=['a,','-,','ab,'] 
        if(cell=="-"):
            output.write("-,")
        else:
            output.write(values[int(cell)])
    if(rowtype=="D2.15b"):
        values=['-,','a,','ab,']
        if(cell=="-"):
            output.write("-,")
        else:
            output.write(values[int(cell)])
    if(rowtype=="B3.7"):
        values=['a,','b,','ab,'] 
        if(cell=="-"):
            output.write("-,")
        else:
            output.write(values[int(cell)])

#het x two null 
    if(rowtype=="D1.11"):
        values=['a,','b,','-,']
        if(cell=="-"):
            output.write("-,") 
        else:
            output.write(values[int(cell)])
    if(rowtype=="D2.16"):
        values=['a,','b,','-,']
        if(cell=="-"):
            output.write("-,")
        else:
            output.write(values[int(cell)])

#het x one null
    if(rowtype=="B1.5"):
        values=['a,','b,','ab,']
        if(cell=="-"):
            output.write("-,") 
        else:
            output.write(values[int(cell)])
    if(rowtype=="B2.6"):
        values=['a,','b,','ab,']
        if(cell=="-"):
            output.write("-,")
        else:
            output.write(values[int(cell)])

# hom x one null
    if(rowtype=="D2.17"):
        values=['a,','-,','ab,']
        if(cell=="-"):
            output.write("-,") 
        else:
            output.write(values[int(cell)])
    if(rowtype=="D2.17b"):
        values=['-,','a,','ab,']
        if(cell=="-"):
            output.write("-,") 
        else:
            output.write(values[int(cell)])
    if(rowtype=="D1.12"):
        values=['a,','-,','ab,']
        if(cell=="-"):
            output.write("-,")
        else:
            output.write(values[int(cell)])
    if(rowtype=="D1.12b"):
        values=['-,','a,','ab,']
        if(cell=="-"):
            output.write("-,")
        else:
            output.write(values[int(cell)])

# accepted missing
    if(rowtype=="A4"):                  ###### LAST FILTERING STEP: separate real A4 from D1.12 / D2.17 if only one of the parents have null alleles
        values=['a,','b,','ab,']
        if(cell=="-"):
            output.write("o,")
        else:
            output.write(values[int(cell)])
    if(rowtype=="A4b"):
        values=['a,','b,','ab,']
        if(cell=="-"):
            output.write("o,")
        else:
            output.write(values[int(cell)])
    if(rowtype=="D1.13"):
        values=['a,','b,','-,']
        if(cell=="-"):
            output.write("o,") 
        else:
            output.write(values[int(cell)])
    if(rowtype=="D2.18"):
        values=['a,','b,','-,']
        if(cell=="-"):
            output.write("o,") 
        else:
            output.write(values[int(cell)])
    if(rowtype=="D1.13x"):
        values=['a,','b,','-,']
        if(cell=="-"):
            output.write("o,") 
        else:
            output.write(values[int(cell)])
    if(rowtype=="D2.18x"):
        values=['a,','b,','-,']
        if(cell=="-"):
            output.write("o,") 
        else:
            output.write(values[int(cell)])

with open('Report_DH21-5921_Altered.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    output=open("26.dec.output.csv","w")
    for row in csv_reader:
        if line_count == 0:
            output.write("filter_PIL,filter_316,filter_CAS,filter3_PIL,filter3_316,filter3_CAS,") 
            output.write("type")
            for cell in row:
                output.write(str(cell)+',')
            output.write('\n')
            line_count += 1
        else:
            output.write(str(check_all_good(row,1))+",")
            output.write(str(check_all_good(row,5))+",")
            output.write(str(check_all_good(row,9))+",")

            typeA=check_3_good(row,1)
            typeB=check_3_good(row,5)
            output.write(str(typeA)+',')
            output.write(str(typeB)+',')
            output.write(str(check_3_good(row,9))+',')
            
            rowtype=get_segtype(typeA,typeB)
            output.write(str(rowtype)+',')

            cell_count=0
            for cell in row:
                if(cell_count<1):
                    output.write(str(cell)+',')
                else:
                    convert_segformat_tooutput(cell,rowtype,output)
                cell_count+=1
            output.write('\n')
            line_count += 1
    output.write("\n"+str(line_count))
    print(str(line_count))
    output.close()
