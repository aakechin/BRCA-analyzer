# This script adds position of mutation in Uljana's reference
# As an input it uses the following arguments:
# (1)inputFile - file with variations
# (2)patients_table - file with table of patients with format
# as patients_table.csv

# v4 - add processing of several mutations at two strands at one position
# v5 - added argparse
# v6 - process of joined avinput file

import glob
import sys
import argparse

def showPercWork(done,allWork):
    percDoneWork=round((done/allWork)*100,1)
    sys.stdout.write("\r"+str(percDoneWork)+"%")
    sys.stdout.flush()

# Readign arguments
par=argparse.ArgumentParser(description="This script adds position of mutation in Uljana's reference")
par.add_argument('--inputFile','-i',dest='inputFile',type=str,help='file with variations',required=True)
par.add_argument('--patientsTable','-p',dest='patientsTable',type=str,help='table with information about each patient: ngs_num patient_id barcode1 barcode2',required=True)
par.add_argument('--torrent','-t',dest='ion_torrent',action='store_true',help='this parameter is needed if you join variants that were called by Torrent Suite Variant Caller')
args=par.parse_args()

file=open(args.inputFile)
pFile=open(args.patientsTable)
resultFile=open(args.inputFile[:-4]+'.withOurCoordinates.xls','w')
ion_torrent=args.ion_torrent
pats={}
for string in pFile:
    cols=string[:-1].split('\t')
    pats[cols[0]]=[cols[1],cols[2]+'_'+cols[3]]

pFile.close()
for string in file:
    if 'Chr\tStart' in string:
        newCols=['Patient_Num','Patient_ID','Barcodes','Chrom','Position',
                 'Position_in_old_ref','Ref','Alt','Qual','Mutation_Type','RefGene_ANN',
                 'Cov_Ref','Cov_Alt','Alt/Total','1000Genomes_All','ExAc_All',
                 'Esp6500_All','Kaviar_AF','dbSNP147','ClinVar_Sign','Cosmic70',
                 'SIFT_Pred','PolyPhen_HDIV_Pred','PolyPhen_HVAR_Pred',
                 'LRT_Pred','MutationTaster_Pred','MutationAccessor_Pred',
                 'FATHMM_Pred','RadialSVM_Pred','LR_Pred']
        resultFile.write('\t'.join(newCols)+'\n')
        continue
    newCols=[]
    cols=string[:-1].split('\t')
    if 'BRCA' not in cols[6]: continue
    cols[0]=cols[0].replace('chr','')
    if cols[0]=='13':
        inURef=int(cols[1])-32889617+1
    elif cols[0]=='17':
        inURef=41277500-int(cols[1])+1
    patients=cols[55].split('|')
    covs=cols[57].split('|')
    patIDs=[]
    barcodes=[]
    refCovs=[]
    altCovs=[]
    altTotals=[]
    for patient,cov in zip(patients,covs):
        if 'patient' in patient:
            if len(patient.split('_'))==2:
                pNum=patient.split('_')[1]
            elif len(patient.split('_'))==3:
                pNum=patient.split('_')[1]+'_'+patient.split('_')[2]
        else:
            pNum=patient
        if pNum in pats.keys():
            patIDs.append(pats[pNum][0])
            barcodes.append(pats[pNum][1])
        else:
            patIDs.append('empty_'+pNum)
            barcodes.append('')
        refCov,altCov=cov.split(',')
        altTotal=str(round(int(altCov)/(int(refCov)+int(altCov)),2))
        refCovs.append(refCov)
        altCovs.append(altCov)
        altTotals.append(altTotal)
    transcriptCheck=False
    if cols[0].replace('chr','')=='13':
        anns=cols[58].split('|')
    else:
        anns=cols[58].split(',')
        for ann in anns:
            a=ann.split('|')
            try:
                if a[6]=='NM_007294.3':
                    anns=a[:]
                    transcriptCheck=True
                    break
            except IndexError as e:
                print('ERROR:',e)
                print(a)
                print(cols)
                exit(0)
        if not transcriptCheck:
            print('ERROR: transcript NM_007294.3 was not found for')
            print(cols)
            exit(0)
    newCols=[cols[55],'|'.join(patIDs),'|'.join(barcodes),cols[0],cols[59],
             str(inURef),cols[60],cols[61],cols[56],anns[1],':'.join([anns[3],anns[6],anns[8],anns[9],anns[10]]),
             '|'.join(refCovs),'|'.join(altCovs),'|'.join(altTotals),cols[23],cols[12],
             cols[11],cols[20],cols[24],cols[25],cols[10],
             cols[31],cols[33],cols[35],
             cols[37],cols[39],cols[41],
             cols[43],cols[45],cols[47]]
    resultFile.write('\t'.join(newCols)+'\n')
file.close()
resultFile.close()
    
