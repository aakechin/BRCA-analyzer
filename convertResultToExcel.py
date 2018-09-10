# This script adds annotation by PredictSNPs
# which uses annotations of SIFT, PolyPhen, MAPP, PhD-SNP and SNAP
# As an input it uses:
# inputFile - file with variations
# siftOutput - file with SIFT output
# pphReportFile - file with PolyPhen output
# annDbDir - directory with annotations files
# outFile - file for output

# Section of importing modules
import sys
import argparse
from multiprocessing import Pool,Queue
import subprocess as sp
import re
import xlsxwriter as xls
import os

# Sections of parameters
QUAL_FILTER2=500
QUAL_FILTER1=1200
ALT_PERC_FILTER1=0.14
ALT_PERC_FILTER2=0.14

# Section of functions
def showPercWork(done,allWork):
    percDoneWork=round((done/allWork)*100,1)
    sys.stdout.write("\r"+str(percDoneWork)+"%")
    sys.stdout.flush()

def convertAminoAcid(aa1):
    aas={'Phe':'F','Leu':'L','Ile':'I','Met':'M','Val':'V','Ser':'S','Pro':'P',
         'Thr':'T','Ala':'A','Tyr':'Y','His':'H','Gln':'Q','Asn':'N','Lys':'K',
         'Asp':'D','Glu':'E','Cys':'C','Trp':'W','Arg':'R','Gly':'G'}
    if aa1 not in aas.keys():
        print('ERROR: No such aminoacid in conversion table!')
        print(aa1)
        exit(0)
    return(aas[aa1])
    
# Section of reading arguments
par=argparse.ArgumentParser(description='This script coverts output to Excel-format')
par.add_argument('--inputFile','-i',dest='inputFile',type=str,help='file with variations')
par.add_argument('--outFile','-o',dest='outFile',type=str,help='file for output')
args=par.parse_args()

thisDir=os.path.dirname(os.path.realpath(__file__))+'/'
wb=xls.Workbook(args.outFile)

### Read file with coordinates for Sanger primers
##coordsFile=open(thisDir+'annotation_databases/brca_sanger_sequencing_primers_coords.xls')
##primersCoords={'17':{},'13':{}}
##for string in coordsFile:
##    cols=string[:-1].split('\t')
##    primersCoords[cols[1]][cols[0]]=[int(cols[2]),int(cols[3])]

# List of SNPs that were already written
file=open(args.inputFile)
checked=[]
wsAll=wb.add_worksheet('All')
wsFiltered=wb.add_worksheet('All_filtered')
wsUnknown=wb.add_worksheet('Unknown')
wsPathogenic=wb.add_worksheet('Pathogenic')
wsPathogenicFiltered=wb.add_worksheet('Pathogenic_filtered')
wsPredictedPathogenic=wb.add_worksheet('Predicted_pathogenic')
# We will filter found variants by the following parameters:
# Alt/total >= 0.14
# For frameshift QUAL >= 2000; for substitutions - QUAL >= 500
f0=wb.add_format({'font_name':'Times New Roman','bold':True,'font_size':12,
                  'align':'center'})
f1=wb.add_format({'font_name':'Times New Roman','bold':False,'font_size':12,
                  'align':'left'})
f2=wb.add_format({'font_name':'Times New Roman','bold':False,'font_size':12,
                  'align':'center'})
colsWidth=[12,15,10,6,10,15,11,11,10,33,
           9,15,18,25,15,15,15,14,10,25,
           14,23,45,26,18,12,22,18,18,21,
           20,21,21,21,21,17,16,18,17]
brca1Domains={'RING':[1,109],'NES':[81,99],'cMyc1':[175,303],'RB':[304,394],'Rad50':[341,748],'cMyc2':[433,511],'NLS1':[503,508],'NLS2':[606,615],
              'Rad51':[758,1064],'Phosph1-CHK2':[988,988],'Phosph2':[1189,1189],'SCD':[1280,1524],'CCD':[1364,1437],'PALB2':[1400,1411],
              'Phosph3-ATM':[1457,1457],'Phosph4-ATM':[1524,1524],'Phosph5-ATM':[1542,1542],'BRCT1':[1650,1745],'BRCT2':[1768,1863]}
brca2Domains={'PALB2':[1,40],'PCAF':[290,453],'NPM1':[639,1000],'BRC1':[1002,1035],'BRC2':[1212,1245],'BRC3':[1421,1454],'BRC4':[1517,1550],
              'BRC5':[1664,1696],'BRC6':[1837,1870],'BRC7':[1972,2004],'BRC8':[2051,2084],'POLH':[1338,1781],'HMG20b':[1648,2190],'FANCD2':[2350,2545],
              'helix':[2479,2667],'DSS1':[2481,2832],'NES':[2682,2698],'OB1':[2682,2794],'OB2':[2804,3054],'OB3':[3073,3167],'NLS1':[3263,3269],
              'NLS2':[3381,3385],'Phosph1':[755,755],'Phosph2':[3291,3291],'Phosph3':[3387,3387],'BubR1':[3189,3418]}
clinVarSignToNums={'Benign':1,'Likely_benign':2,'Conflicting_interpretations_of_pathogenicity':3,'Uncertain_significance':3,'other':3,'.':3,'not_provided':3,'Likely_pathogenic':4,'Pathogenic':5}
i1=0; i2=0; i3=0; i4=0; i5=0; i6=0
pp=re.compile('\d+')
for string in file:
    cols=string[:-1].split('\t')
    if 'Chrom\tPosition' in string:
        newCols=cols[:10]+['Gene_Name','Transcript_ID','Exon_Num','CDS_HGVS','Protein_HGVS','Domain']+cols[11:]
        for k,newCol in enumerate(newCols):
            # formula for width was calculated by regression
            wsAll.set_column(k,k,colsWidth[k])
            wsFiltered.set_column(k,k,colsWidth[k])
            wsUnknown.set_column(k,k,colsWidth[k])
            wsPathogenic.set_column(k,k,colsWidth[k])
            wsPathogenicFiltered.set_column(k,k,colsWidth[k])
            wsPredictedPathogenic.set_column(k,k,colsWidth[k])
        wsAll.write_row(0,0,newCols,f0)
        wsFiltered.write_row(0,0,newCols,f0)
        wsUnknown.write_row(0,0,newCols,f0)
        wsPathogenic.write_row(0,0,newCols,f0)
        wsPathogenicFiltered.write_row(0,0,newCols,f0)
        wsPredictedPathogenic.write_row(0,0,newCols,f0)
        i1+=1; i2+=1; i3+=1; i4+=1; i5+=1; i6+=1
        continue
    transcriptCheck=False
    # Annotate mutations that influence on the protein sequence by domain
    domains=[]
    if cols[10]!='0':
        ann=cols[10].split(':')
        if ann[4]!='':
            try:
                pm=pp.findall(ann[4])
            except IndexError as e:
                print('ERROR:',e)
                print(ann)
                exit(0)
            if len(pm)==0:
                print('ERROR: No found numbers in Protein designation!',ann[4])
                print(ann)
                exit(0)
            if ann[0]=='BRCA1':
                for key,value in brca1Domains.items():
                    if int(pm[0])>=value[0] and int(pm[0])<=value[1]:
                        domains.append(key)
            elif ann[0]=='BRCA2':
                for key,value in brca2Domains.items():
                    if int(pm[0])>=value[0] and int(pm[0])<=value[1]:
                        domains.append(key)
            else:
                print('ERROR: unknown gene!',anns[3])
                exit(0)
            newCols=cols[:10]+ann+['&'.join(domains)]+cols[11:]
        else:
            newCols=cols[:10]+ann+['']+cols[11:]
    else:
        newCols=cols[:10]+['','','','','','']+cols[11:]
    patNums=cols[0].split('|')
    patIds=cols[1].split('|')
    barcodes=cols[2].split('|')
    quals=cols[8].split('|')
    covRefs=cols[11].split('|')
    covAlts=cols[12].split('|')
    altTotals=cols[13].split('|')
    clinVarSigns=newCols[27].split('|')
    clinVarSum=[]
    for clin in clinVarSigns:
        clins=clin.split('/')
        for cl in clins:
            try:
                clinVarSum.append(clinVarSignToNums[cl])
            except KeyError:
                print('ERROR!',clin)
                print(string)
                exit(0)
    clinVarSign=int(round(sum(clinVarSum)/len(clinVarSum),0))
    for patNum,patId,barcode,qual,covRef,covAlt,altTotal in zip(patNums,patIds,barcodes,quals,covRefs,covAlts,altTotals):
        newCols[0]=patNum; newCols[1]=patId; newCols[2]=barcode; newCols[8]=qual
        newCols[16]=covRef; newCols[17]=covAlt; newCols[18]=altTotal
        # Write all variants to worksheet "All"
        wsAll.write_row(i1,0,newCols[0:2],f1)
        wsAll.write_row(i1,2,newCols[2:],f2)
        i1+=1
    ##    print(newCols[9],newCols[8],newCols[18],newCols[19],newCols[21],newCols[24])
        # Write all filtered
        if float(newCols[18])>=ALT_PERC_FILTER1 or ((float(newCols[18])>=ALT_PERC_FILTER2) and float(newCols[8])>=QUAL_FILTER2) or float(newCols[8])>=QUAL_FILTER1:
            wsFiltered.write_row(i2,0,newCols[0:2],f1)
            wsFiltered.write_row(i2,2,newCols[2:],f2)
            i2+=1
        # Write pathogenic
        if ((('frameshift_variant' in newCols[9] or 'stop_gained' in newCols[9] or 'stop_lost' in newCols[9]
            or 'start_lost' in newCols[9] or 'splice_acceptor_variant' in newCols[9]
            or 'splice_donor_variant' in newCols[9]) and newCols[25]!='no' and clinVarSign>=3) or newCols[25]=='yes' or clinVarSign>=4):
            wsPathogenic.write_row(i3,0,newCols[0:2],f1)
            wsPathogenic.write_row(i3,2,newCols[2:],f2)
            i3+=1
            # Write pathogenic filtered
            if float(newCols[18])>=ALT_PERC_FILTER1 or ((float(newCols[18])>=ALT_PERC_FILTER2) and float(newCols[8])>=QUAL_FILTER2) or float(newCols[8])>=QUAL_FILTER1:
                wsPathogenicFiltered.write_row(i4,0,newCols[0:2],f1)
                wsPathogenicFiltered.write_row(i4,2,newCols[2:],f2)
                i4+=1
        # Write predicted pathogenic filtered
        # Low frequency in 1000Genomes, ExAc, ESP6500 and Kaviar, it's not Benign by ClinVar and BIC
        # and it's predicted by all tools as deleterious
        elif (cols[24]=='D' and (cols[25]=='D' or cols[25]=='P') and (cols[26]=='D' or cols[26]=='P') and cols[27]=='D' and (cols[28]=='A' or cols[28]=='D')
              and (cols[29]=='H' or cols[29]=='M') and cols[30:]==['D','D','D']
              and (cols[20]=='' or cols[20]=='unknown') and clinVarSign>=3):
            wsPredictedPathogenic.write_row(i5,0,newCols[0:2],f1)
            wsPredictedPathogenic.write_row(i5,2,newCols[2:],f2)
            i5+=1
        # Write unknown
        elif ((newCols[19]=='.' or float(newCols[19])<0.005) and newCols[25]!='no' and clinVarSign==3 and len(patNums)<5
              and newCols[9]!='synonymous_variant' and newCols[9]!='downstream_gene_variant' and newCols[9]!='upstream_gene_variant'):
            # If it is an intron variant, then we check if it is not more than 50 nucleotides before 3'-end of intron
            ## 50 nucleotidesm of 3'-end of intron is a possible polypirimidine tract
            intronVarPat=re.compile('c.\d+\-(\d+)')
            if newCols[9]=='intron_variant':
                intronVarMatches=intronVarPat.findall(newCols[13])
                if len(intronVarMatches)>0 and int(intronVarMatches[0])<=50:
                    wsUnknown.write_row(i6,0,newCols[0:2],f1)
                    wsUnknown.write_row(i6,2,newCols[2:],f2)
                    i6+=1
            else:
                wsUnknown.write_row(i6,0,newCols[0:2],f1)
                wsUnknown.write_row(i6,2,newCols[2:],f2)
                i6+=1
wb.close()
file.close()
