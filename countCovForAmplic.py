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


import os
import glob
import sys
import subprocess as sp
import re
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
    output=sp.check_output(configs[1]+"samtools mpileup -Q 0 -d 1000000 -r chr17:41196312-41277500 -t DP -A "+bamFile+" > "+resultFile+'2',shell=True,stderr=sp.STDOUT)
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

# Make file with coordinates
# Read file with primers
minMean=sys.argv[1]
coordsFileName=sys.argv[2]
patFileName=sys.argv[3]
bamFilesSpec=sys.argv[4]
resultFileName=sys.argv[5]
threads=sys.argv[6]

if minMean not in ['min','mean','max']:
    print('ERROR! The last argument should be "min", "mean" or "max"')
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
pFile=open(sys.argv[3])
pats={}
for string in pFile:
    cols=string[:-1].split('\t')
    pats[cols[0]]=[cols[1],cols[2]+'_'+cols[3]]
pFile.close()

resultFile=open(resultFileName,'w')
resultFile.write("Patient#\tPatient_ID\tBarcodes\tMedian_Coverage\tNumber_<30")
for k,exon in enumerate(coords):
    resultFile.write("\tamplicon#%d" % (k+1))
resultFile.write("\n")
allWork=len(bamFiles)
showPercWork(0,allWork)
p=Pool(threads)
doneWork=0
for res in p.imap_unordered(processBamFile,bamFiles,10):
    bamFile,covs=res
    bamFilePart=bamFile[bamFile.rfind('/')+1:]
    resultFile.write("%s" % (bamFilePart[:bamFilePart.index('.')].replace('patient_','')))
    resultFile.write('\t'+'\t'.join(pats[bamFilePart[:bamFilePart.index('.')].replace('patient_','')]))
    resultFile.write('\t'+'\t'.join(covs)+'\n')
    doneWork+=1
    showPercWork(doneWork,allWork)
##for i,bamFile in enumerate(bamFiles):
##    mpileupFile=bamFile+".mpileup"
##    mpileup(bamFile,mpileupFile)
##    bamFilePart=bamFile[bamFile.rfind('/')+1:]
##    resultFile.write("%s" % (bamFilePart[:bamFilePart.index('.')].replace('patient_','')))
##    resultFile.write('\t'+'\t'.join(pats[bamFilePart[:bamFilePart.index('.')].replace('patient_','')]))
##    covs=[]
##    for k,coord in enumerate(coords):
##        if coord[0]==13:
##            cov=round(countCoverage(bamFile,mpileupFile+'1',coord[0],coord[1],coord[2],minMean),3)
##        elif coord[0]==17:
##            cov=round(countCoverage(bamFile,mpileupFile+'2',coord[0],coord[1],coord[2],minMean),3)
##        covs.append(cov)
##    covLess30=sum(z<30 for z in covs)
##    covs=list(map(str,[round(median(covs),3),covLess30]+covs))
##    resultFile.write('\t'+'\t'.join(covs)+'\n')
##    showPercWork(i+1,allWork)
resultFile.close()
print()
