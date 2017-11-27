# This script gets percent of properly trimmed reads

import glob,re,argparse,sys
import subprocess as sp

def showPercWork(done,allWork):
    percDoneWork=round((done/allWork)*100,2)
    sys.stdout.write("\r"+str(percDoneWork)+"%")
    sys.stdout.flush()

par=argparse.ArgumentParser(description='This script gets percent of properly trimmed reads')
par.add_argument('--native-reads','-nat',dest='nativeReadsFiles',type=str,help='regular expression for files with native R1 reads (better with Undertermined reads)',required=True)
par.add_argument('--trimmed-reads','-trim',dest='trimReadsFiles',type=str,help='regular expression for files with trimmed R1 reads. Use this argument, only if have trimmed reads',required=False)
par.add_argument('--patients-file','-pat',dest='patFile',type=str,help='file with patients information (Patient_Num, Patient_ID, index1, index2)',required=True)
par.add_argument('--output-file','-out',dest='outFile',type=str,help='file for output',required=True)
args=par.parse_args()

file=open(args.patFile)
pats={}
for string in file:
    if 'PatientID' in string: continue
    cols=string.replace('\n','').split('\t')
    pats[cols[0]]=cols[1:4]

p1=re.compile('(\d+)_S\d+_')
p2=re.compile('patient_(\d+).r1')
nump=re.compile('(\d+)')
ds1=glob.glob(args.nativeReadsFiles)
reads1Num={}
if args.trimReadsFiles:
    ds2=glob.glob(args.trimReadsFiles)
    reads2Num={}
allWork=len(ds1)
print('Evaluating number of reads for each sample...')
showPercWork(0,allWork)
for i,d in enumerate(sorted(ds1)):
    if 'Undetermined' in d:
        patNum='Undetermined'
    else:
        patNum=p1.findall(d)[0]
    out=sp.check_output('zcat '+d.replace(' ','\ ')+' | wc -l',shell=True,stderr=sp.STDOUT,universal_newlines=True)
    totalReadsNum=int(nump.findall(out)[0])
    reads1Num[patNum]=totalReadsNum
    showPercWork(i+1,allWork)
print()
if args.trimReadsFiles:
    allWork=len(ds2)
    print('Evaluating number of properly trimmed reads...')
    showPercWork(0,allWork)
    for i,d in enumerate(sorted(ds2)):
        patNum=p2.findall(d)[0]
        out=sp.check_output('zcat '+d.replace(' ','\ ')+' | wc -l',shell=True,stderr=sp.STDOUT,universal_newlines=True)
        trimmedReadsNum=int(nump.findall(out)[0])
        reads2Num[patNum]=trimmedReadsNum
        showPercWork(i+1,allWork)
    print()
rFile=open(args.outFile,'w')
if args.trimReadsFiles:
    rFile.write('\t'.join(['Patient_Num','Patient_ID','Barcode_1','Barcode_2','Total_Reads','Properly_Trimmed','Share_of_Reads'])+'\n')
    for key,item in reads1Num.items():
        if key not in reads2Num and key!='Undetermined':
            perc=0; reads2Num[key]=0
        elif key=='Undetermined':
            rFile.write('\t'.join([key,'','','',str(reads1Num[key]),'',''])+'\n')
            continue
        else:
            perc=reads2Num[key]/reads1Num[key]
        rFile.write('\t'.join([key,pats[key][0],pats[key][1],pats[key][2],str(reads1Num[key]),str(reads2Num[key]),str(perc)])+'\n')
else:
    rFile.write('\t'.join(['Patient_Num','Patient_ID','Barcode_1','Barcode_2','Total_Reads'])+'\n')
    for key,item in reads1Num.items():
        if key=='Undetermined':
            rFile.write('\t'.join([key,'','','',str(reads1Num[key])])+'\n')
            continue
        rFile.write('\t'.join([key,pats[key][0],pats[key][1],pats[key][2],str(reads1Num[key])])+'\n')
rFile.close()