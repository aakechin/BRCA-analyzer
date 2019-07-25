# VERY IMPORTANT!
# mpileup of samtools count only reads that are paired
# at the same time UGENE count all reads including ones that were not paired
# It uses the followng input arguments:
# [min or mean] - which coverage to output: minimal for amplicon bases or an average one
# [coordsFileName] - file for storing coordinates of each primer
# [patientsTable] - file with table of patients information
# [bamFilesSpec] - regular expression of BAM-files
# [resultFileName] - file for results
# [threads]


import os,glob,sys,re,argparse
import subprocess as sp
from statistics import mean
from statistics import median
from multiprocessing import Pool

global thisDir,configs,coords,minMean
thisDir=os.path.dirname(os.path.realpath(__file__))+'/'
configs=open(thisDir+'config.txt').read().split('\n')

def showPercWork(done,allWork):
    percDoneWork=round((done/allWork)*100,1)
    sys.stdout.write("\r"+str(percDoneWork)+"%")
    sys.stdout.flush()

# Next function counts average read depth for some region of particular chromosome
def countCoverage(bamFile,mpileupFile,chrom=17,start=0,end=0,minMean='mean'):
    if end<start:
        t=start
        start=end
        end=t
    file=open(mpileupFile)
    covs=[] # coverages
    prevCoord=start-1
    for l in file:
        cols=l[:-1].split('\t')
        coord=int(cols[1])
        if cols[0].replace('chr','')=='X': continue
        if chrom!=int(cols[0].replace('chr','')) or coord<start:
            continue
        elif end>0 and int(cols[1])>end:
            break
        else:
            if coord>prevCoord+1:
                covs.extend([0]*(coord-prevCoord-1))
            covs.append(int(cols[3]))
            prevCoord=coord
    # Get last position
    if len(covs)==0:
        return 0
    else:
        if minMean=='mean': return mean(covs)
        elif minMean=='min': return min(covs)
        elif minMean=='max': return max(covs)
    
# Next function creates file mpileup with coverages for each position of BAM-file
def mpileup(bamFile,resultFile="output.mpileup"):
    output=sp.check_output(configs[1]+"samtools mpileup -Q 0 -d 1000000 -r chr13:32889617-32973809 -t DP -A "+bamFile+" > "+resultFile+'1',shell=True,stderr=sp.STDOUT)
    output=sp.check_output(configs[1]+"samtools mpileup -Q 0 -d 1000000 -r chr17:41196312-41279700 -t DP -A "+bamFile+" > "+resultFile+'2',shell=True,stderr=sp.STDOUT)
    return resultFile

def processBamFile(bamFile):
    mpileupFile=bamFile+".mpileup"
    mpileup(bamFile,mpileupFile)
    covs=[]
    for k,coord in enumerate(coords):
        if coord[0]==13:
            cov=round(countCoverage(bamFile,mpileupFile+'1',coord[0],coord[1],coord[2],minMean),3)
        elif coord[0]==17:
            cov=round(countCoverage(bamFile,mpileupFile+'2',coord[0],coord[1],coord[2],minMean),3)
        covs.append(cov)
    covLess30=sum(z<30 for z in covs)
    covs=list(map(str,[round(median(covs),3),covLess30]+covs))
    return(bamFile,covs)

# Readign arguments
par=argparse.ArgumentParser(description="This script evaluates coverage of amplicons")
par.add_argument('--min-mean-max','-mmm',dest='minMeanMax',type=str,help='min, mean or max values should be outputed by the program for each amplicon. Default: mean',required=False,default='mean')
par.add_argument('--file-with-coordinates','-coord',dest='coordsFileName',type=str,help='TSV-file with coordinates for each amplicon: amplicon_name | chromosome | start | end',required=True)
par.add_argument('--pateints-file','-pat',dest='patFileName',type=str,help='TSV-file with information about each sample: sample_number | sample ID | index1| index2',required=False)
par.add_argument('--bam-files','-bam',dest='bamFilesSpec',type=str,help='regular expression for BAM-files',required=True)
par.add_argument('--result-file-name','-res',dest='resultFileName',type=str,help='name for file with results',required=True)
par.add_argument('--patientsList','-pl',dest='patientsList',type=str,help='list of sample numbers that correspond IDs of files with reads (BRCA-alayzer send it to this tool automatically)',required=False)
par.add_argument('--threads','-th',dest='threads',type=int,help='number of threads',default=2)
args=par.parse_args()

# Make file with coordinates
# Read file with primers
minMean=args.minMeanMax
coordsFileName=args.coordsFileName
patFileName=args.patFileName
bamFilesSpec=args.bamFilesSpec
resultFileName=args.resultFileName
threads=args.threads
if args.patientsList:
    patNums=args.patientsList.split('_')
else:
    patNums=None 

if minMean not in ['min','mean','max']:
    print('ERROR! The argument -mmm (--min-mean-max) should be "min", "mean" or "max"')
    exit(0)
bamFiles=glob.glob(bamFilesSpec)
if len(bamFiles)==0:
    print('ERROR: There is no BAM-files selected')
    print(sys.argv)
    exit(0)
elif len(bamFiles)==1:
    print('There was only one BAM-file selected.')
try:
    threads=int(threads)
except ValueError:
    print('ERROR: Number of threads should be integer. You write:\n',threads);
    exit(0)

# Read file with coordinates
file=open(coordsFileName)
coords=[]
for string in file:
    cols=string[:-1].split('\t')
    coords.append([int(cols[1].replace('chr','')),int(cols[2]),int(cols[3])])

# Read file with patients table
pats={}
if patFileName:
    pFile=open(patFileName)
    for string in pFile:
        if ('PatientID' in string or
            'Patient_Num' in string or
            'Patient_ID' in string or
            string=='' or
            string==' ' or
            string=='\n'):
            continue
        cols=string.replace('\n','').replace('\r','').split('\t')
        pats[cols[0]]=[cols[1],cols[2]+'_'+cols[3]]
    pFile.close()
else:
    for pat in patNums:
        pats[pat]=['N/A','N/A_N/A']

resultFile=open(resultFileName,'w')
resultFile.write("Patient#\tPatient_ID\tBarcodes\tMedian_Coverage\tNumber_<30")
for k,exon in enumerate(coords):
    resultFile.write("\tamplicon#%d" % (k+1))
resultFile.write("\n")
allWork=len(bamFiles)
showPercWork(0,allWork)
p=Pool(threads)
doneWork=0
showPercWork(doneWork,allWork)
for res in p.imap_unordered(processBamFile,bamFiles,10):
    bamFile,covs=res
    bamFilePart=bamFile[bamFile.rfind('/')+1:]
    resultFile.write("%s" % (bamFilePart[:bamFilePart.index('.')].replace('patient_','')))
    patNum=bamFilePart[:bamFilePart.index('.')].replace('patient_','')
    if patNum in pats.keys():
        resultFile.write('\t'+'\t'.join(pats[patNum]))
    elif patNums and patNum in patNums:
        resultFile.write('\t'+'\t'.join(pats[str(patNums.index(patNum)+1)]))
    else:
        resultFile.write('\t'*2)
    resultFile.write('\t'+'\t'.join(covs)+'\n')
    doneWork+=1
    showPercWork(doneWork,allWork)
resultFile.close()
print()
