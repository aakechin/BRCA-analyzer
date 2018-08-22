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
par.add_argument('--patients-file','-pat',dest='patFile',type=str,help='file with patients information (Patient_Num, Patient_ID, index1, index2)',required=False)
par.add_argument('--output-file','-out',dest='outFile',type=str,help='file for output',required=True)
args=par.parse_args()

if args.patFile:
    file=open(args.patFile)
    pats={}
    for string in file:
        if 'PatientID' in string: continue
        cols=string.replace('\n','').split('\t')
        pats[cols[0]]=cols[1:4]

p1=re.compile('(\d+)')
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
filesPatNums=[]
for i,d in enumerate(sorted(ds1)):
    if 'Undetermined' in d:
        patNum='Undetermined'
    else:
        dName=d[d.rfind('/')+1:]
        patNum=p1.findall(dName)[0]
        filesPatNums.append(int(patNum))
    out=sp.check_output('zcat '+d.replace(' ','\ ')+' | wc -l',shell=True,stderr=sp.STDOUT,universal_newlines=True)
    totalReadsNum=int(nump.findall(out)[0])
    reads1Num[patNum]=totalReadsNum/2
    showPercWork(i+1,allWork)
filesPatNums=sorted(filesPatNums)
##print(reads1Num)
print()
if args.trimReadsFiles:
    allWork=len(ds2)
    print('Evaluating number of properly trimmed reads...')
    showPercWork(0,allWork)
    for i,d in enumerate(sorted(ds2)):
        m2=p2.findall(d)
        if len(m2)==0: patNum=p1.findall(d)[0]
        else: patNum=m2[0]
        out=sp.check_output('zcat '+d.replace(' ','\ ')+' | wc -l',shell=True,stderr=sp.STDOUT,universal_newlines=True)
        trimmedReadsNum=int(nump.findall(out)[0])
        reads2Num[patNum]=trimmedReadsNum/2
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
        if args.patFile:
            if key in pats.keys():
                rFile.write('\t'.join([key,pats[key][0],pats[key][1],pats[key][2],str(reads1Num[key]),str(reads2Num[key]),str(perc)])+'\n')
            else:
                rFile.write('\t'.join([key,'empty_'+key,'-','-',str(reads1Num[key]),str(reads2Num[key]),str(perc)])+'\n')
        else:
            rFile.write('\t'.join([key,'N/A','N/A','N/A',str(reads1Num[key]),str(reads2Num[key]),str(perc)])+'\n')
else:
    rFile.write('\t'.join(['Patient_Num','Patient_ID','Barcode_1','Barcode_2','Total_Reads'])+'\n')
    for key,item in reads1Num.items():
        if key=='Undetermined':
            rFile.write('\t'.join([key,'','','',str(reads1Num[key])])+'\n')
            continue
        if args.patFile:
            if key not in pats.keys():
                if str(filesPatNums.index(int(key))+1) in pats.keys():
                    newKey=str(filesPatNums.index(int(key))+1)
                    print("WARNING: Patients' numbers in patients table do not correspond numbers in the file\n"
                          'So the program will take them as an order numbers')
                    rFile.write('\t'.join([key,pats[newKey][0],pats[newKey][1],pats[newKey][2],str(reads1Num[key])])+'\n')
                else:
                    print("WARNING: Patients' numbers in patients table do not correspond numbers in the file\n"
                          'So the program will take this sample without patient ID')
                    rFile.write('\t'.join([key,'','','',str(reads1Num[key])])+'\n')
            else:
                 rFile.write('\t'.join([key,pats[key][0],pats[key][1],pats[key][2],str(reads1Num[key])])+'\n')
        else:
            rFile.write('\t'.join([key,'N/A','N/A','N/A',str(reads1Num[key])])+'\n')
rFile.close()
