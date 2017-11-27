# This script joins mutation of all patients to one file
# Use the following arguments:
# [1] - regular expression for variants files
# [2] - name for resultFile
## v4 - added argparse
## v5 - join mutations of avinput (annovar input)

import glob
import sys
import argparse
import re

def showPercWork(done,allWork):
    import sys
    percDoneWork=round((done/allWork)*100,1)
    sys.stdout.write("\r"+str(percDoneWork)+"%")
    sys.stdout.flush()

# Readign arguments
par=argparse.ArgumentParser(description='This script joins mutation of all patients to one file')
par.add_argument('--varFiles','-v',dest='varFiles',type=str,help='regular expression for choosing files with variations',required=True)
par.add_argument('--outFile','-o',dest='outFile',type=str,help='directory for output',required=True)
args=par.parse_args()

ds=glob.glob(args.varFiles)
if len(ds)==0:
    print('ERROR: No files were selected! Maybe you write "~/" but you should use /home/USERNAME/')
    print(args.varFiles)
    exit(0)
elif len(ds)==1:
    print('WARNING: Only one file was selected!')
    print(args.varFiles)

titlesCheck=False
allWork=len(ds)
showPercWork(0,allWork)
allData={}
annp=re.compile('ANN\=([^\;]+)')
for i,d in enumerate(ds):
    file=open(d)
    dPart=d[d.rfind('/')+1:]
    patName=dPart[:dPart.index('.')]
    adi=1
    for string in file:
        if 'CHROM' in string:
            if not titlesCheck:
                titlesCheck=True
            continue
        cols=string[:-1].split('\t')
        qual=cols[10]
        ads=cols[14].split(':')[1].split(',')
        adNum=len(ads)
        # If we meet locus with several alleles first time
        if adNum>2:
            ad=ads[0]+','+ads[adi]
            if adi<adNum-1:
                adi+=1
            elif adi==adNum-1:
                adi=1
        else:
            ad=','.join(ads)
        annm=annp.findall(cols[12])
        ann=annm[0]
        pos=cols[6]
        ref=cols[8]
        alt=cols[9]
        key='\t'.join(cols[:5])
        if key not in allData.keys():
            allData[key]=[patName,qual,ad,ann,pos,ref,alt]
        else:
            allData[key][0]+='|'+patName
            allData[key][1]+='|'+qual
            allData[key][2]+='|'+ad
    file.close()
    showPercWork(i+1,allWork)
resultFile=open(args.outFile,'w')
for key,item in allData.items():
    resultFile.write(key+'\t'+'\t'.join(item)+'\n')
resultFile.close()
print()
        
