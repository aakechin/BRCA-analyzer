# This script takes annotation from 1000Genomes, ClinVar, dbSNP and BIC
# The result file is like inputFile but with word 'clinical'
# As an input argument it uses:
# (1)inputFile - file with variations
## v3 - added argparse
## v4 - annotates after ANNOVAR. BIC annotation adds trivial BIC designation

# Section of importing modules
import re
from Bio import SeqIO
from Bio import pairwise2
import sys,os
import argparse

# Global variables
global thisDir
thisDir=os.path.dirname(os.path.realpath(__file__))+'/'

# Section of functions
def checkNum(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def showPercWork(done,allWork):
    percDoneWork=round((done/allWork)*100,1)
    sys.stdout.write("\r"+str(percDoneWork)+"%")
    sys.stdout.flush()

def readBicDb(brca):
    # brca - is a string of number - '1' or '2'
    file=open(thisDir+'annotation_databases/brca'+brca+'_BIC_data.txt')
    lines=file.read().split('\n')
    gPosP=re.compile('(\d+)(?:_(\d+))?')
    gPosSubstP=re.compile('([ATGC]+)>([ATGC]+)')
    gPosInDelP=re.compile('(del|ins)([ATGC\d]+)?')
    # We make list of lists that contain the following information:
    # [varType, posInDb, refAllele, altAllele, regionStart(-200), regionEnd(+200), cDnaPos, clinImp, mutAlleleSeq]
    bicDB=[]
    bicCoords=[]
    substCheck=[]
    indelCheck=[]
    for l in lines[:-1]:
        if 'HGVS Genomic (hg19)' in l:
            cols=l.split('\t')
            trivDesignNum=cols.index('Designation')
            cDnaColNum=cols.index('HGVS cDNA')
            gPosColNum=cols.index('HGVS Genomic (hg19)')
            clinImpColNum=cols.index('Clinically Importance')
            continue
        cols=l.split('\t')
        try:
            cDnaPos=cols[cDnaColNum]
        except IndexError as e:
            print('ERROR:',e)
            print(cols)
            print(l)
            exit(0)
        gPos=cols[gPosColNum]
        clinImp=cols[clinImpColNum]
        trivDesign=cols[trivDesignNum]
        gPosM=gPosP.findall(gPos)
        gPosSubstM=gPosSubstP.findall(gPos)
        gPosInDelM=gPosInDelP.findall(gPos)
        if gPos=='-':
            continue
        # If this variation is substitution
        if len(gPosSubstM)>0:
            if [gPosM[0][0],gPosSubstM[0][0],gPosSubstM[0][1]] not in substCheck:
                bicDB.append(['substitution',gPosM[0][0],gPosSubstM[0][0],gPosSubstM[0][1],0,0,cDnaPos,trivDesign,clinImp,''])
                bicCoords.append(int(gPosM[0][0]))
                substCheck.append([gPosM[0][0],gPosSubstM[0][0],gPosSubstM[0][1]])
                continue
            else: continue
        # If this is not simple substitution, we create mutated sequence of this allele
        try:
            if gPosM[0][1]=='':
                regionStart=int(gPosM[0][0])-200
                regionEnd=int(gPosM[0][0])+200
            else:
                regionStart=int(gPosM[0][0])-200
                regionEnd=int(gPosM[0][1])+200
        except IndexError as e:
            print('ERROR:',e)
            print(gPosM)
            print(gPos)
            exit(0)
        # If there is insertions and/or delitions in variation
        # There are only the following variants of InDels:
        # del
        # ins
        # del ins
        # If insertion is the first
        if len(gPosInDelM)==1:
            indel=gPosInDelM[0]
            # The first - only ins
            if gPosInDelM[0][0]=='ins':
                altAlleleSeq=refSeqs[brca][regionStart-1:int(gPosM[0][0])]+indel[1]+refSeqs[brca][int(gPosM[0][0]):regionEnd]
                if altAlleleSeq not in indelCheck:
                    bicDB.append(['insertion',gPosM[0][0],'-',indel[1],regionStart,regionEnd,cDnaPos,trivDesign,clinImp,altAlleleSeq])
                    bicCoords.append(int(gPosM[0][0]))
                    indelCheck.append(altAlleleSeq)
            # The second - only del
            elif gPosInDelM[0][0]=='del':
                # This variant is because of idiots who designate deletion with numbers: del213!
                if indel[1]=='' or not ('A' in indel[1] or 'T' in indel[1] or 'C' in indel[1] or 'G' in indel[1]):
                    altAlleleSeq=refSeqs[brca][regionStart-1:int(gPosM[0][0])-1]+refSeqs[brca][int(gPosM[0][1]):regionEnd]
                else:
                    altAlleleSeq=refSeqs[brca][regionStart-1:int(gPosM[0][0])-1]+refSeqs[brca][int(gPosM[0][0])+len(indel[1])-1:regionEnd]
                if altAlleleSeq not in indelCheck:
                    bicDB.append(['deletion',gPosM[0][0],indel[1],'-',regionStart,regionEnd,cDnaPos,trivDesign,clinImp,altAlleleSeq])
                    bicCoords.append(int(gPosM[0][0]))
                    indelCheck.append(altAlleleSeq)
        # The third del+ins
        # Additional check, that there is not another variants
        elif len(gPosInDelM)==2 and gPosInDelM[0][0]=='del' and gPosInDelM[1][0]=='ins':
            if altAlleleSeq not in indelCheck:
                indel=gPosInDelM[0]
                altAlleleSeq=refSeqs[brca][:int(gPosM[0][0])-1]+refSeqs[brca][int(gPosM[0][0])-1+len(indel[1]):]
                indel=gPosInDelM[1]
                altAlleleSeq=altAlleleSeq[regionStart-1:int(gPosM[0][0])-1]+indel[1]+altAlleleSeq[int(gPosM[0][0])-1:regionEnd]
                bicCoords.append(int(gPosM[0][0]))
                bicDB.append(['deletion-insertion',gPosM[0][0],gPosInDelM[0][1],gPosInDelM[1][1],regionStart,regionEnd,cDnaPos,trivDesign,clinImp,altAlleleSeq])
                indelCheck.append(altAlleleSeq)
        else:
            print('INTERNAL ERROR: there is no such type of InDel')
            print(gPos)
            exit(0)
    bicDB.sort(key=lambda x:int(x[1]))
    bicCoords.sort()
    return [bicDB,bicCoords]

# Reading arguments
par=argparse.ArgumentParser(description='This script takes annotation from 1000Genomes, ClinVar, dbSNP and BIC')
par.add_argument('--inputFile','-i',dest='inputFile',type=str,help='file with variations',required=True)
args=par.parse_args()

# Read fasta sequences of chr13 and chr17
global refSeqs
refSeqs={}
refData=SeqIO.parse(thisDir+'ref/human_g1k_v37_chr13+17.fasta','fasta')
# refSeqs[0] - chr13, refSeqs[1] - chr17
k=2
for f in refData:
    refSeqs[str(k)]=str(f.seq)
    k-=1

brca_bic_db={}
brca_bic_db['17']=readBicDb('1')
brca_bic_db['13']=readBicDb('2')
chrToBrca={'13':'2','17':'1'}

inputFile=open(args.inputFile)
resultFile=open(args.inputFile[:-4]+'.clinical.xls','w')
lines=inputFile.read().split('\n')
inputFile.close()
allWork=len(lines[:-1])
showPercWork(0,allWork)
for i,l in enumerate(lines[:-1]):
    if 'Chrom\tPosition' in l:
        cols=l.split('\t')
        newCols=cols[:18]+['CDS_BIC','Trivial_BIC','BIC_Sign']+cols[18:]
        resultFile.write('\t'.join(newCols)+'\n')
        continue
    cols=l.split('\t')
    newCols=[]
    pos=int(cols[4])
    refAl=cols[6]
    altAl=cols[7]
    brca=chrToBrca[cols[3]]
    # Annotation by BIC
    bicCheck=False
    # First of all we determine type of variation
    # if it is a substitution, we try to find it by coordinate
    if refAl!='-' and altAl!='-':
        # we find it in the list of coordinates
        if pos in brca_bic_db[cols[3]][1]:
            posStart=brca_bic_db[cols[3]][1].index(pos)
            for b in brca_bic_db[cols[3]][0][posStart:]:
                if refAl==b[2] and altAl==b[3]:
                    bicCds=b[6]
                    bicClin=b[7]
                    bicDesign=b[8]
                    newCols=cols[:18]+b[6:9]+cols[18:]
                    break
                # If we came to the line where positions are not equal, break it
                elif str(pos)!=b[1]:
                    bicCds='-'
                    bicClin='-'
                    bicDesign='-'
                    newCols=cols[:18]+['','','']+cols[18:]
                    break
        # If we didn't manage to find position in list of coordinates
        else:
            bicCds='-'
            bicClin='-'
            bicDesign='-'
            newCols=cols[:18]+['','','']+cols[18:]
    # if it is insertion
    elif refAl=='-':
        # we convert it to the BIC format, do not increasing value of position
        # because now we annotates with ANNOVAR
        # and make sequence
        insNucs=altAl
        # make mutated sequence region
        regionStart=pos-200
        regionEnd=pos+200
        altAlleleSeq=refSeqs[brca][regionStart-1:pos]+insNucs+refSeqs[brca][pos:regionEnd]
        # we search coordinate from BIC that is not less far than 200 bp
        for b in brca_bic_db[cols[3]][0]:
            # if we found such coordinate
            if pos-int(b[1])<=50 and b[0]=='insertion':
                # we align our mutated sequence with seq from db
                a=pairwise2.align.globalms(b[8],altAlleleSeq,2,-1,-1.53,0)
                # This alignment has the following format:
                # [('---TCAGTAACAAATGCTCCTATA-', 'GGCTCAGTAACAAATGCTCCTATAA', 21.0, 0, 25),...]
                # we calculate offset for both sides of both sequences
                ls1=len(a[0][0])-len(a[0][0].lstrip('-'))
                ls2=len(a[0][1])-len(a[0][1].lstrip('-'))
                rs1=len(a[0][0].rstrip('-'))
                rs2=len(a[0][1].rstrip('-'))
                if a[0][0][max(ls1,ls2):min(rs1,rs2)]==a[0][1][max(ls1,ls2):min(rs1,rs2)]:
                    bicCds=b[6]
                    bicClin=b[7]
                    bicDesign=b[8]
                    newCols=cols[:18]+b[6:9]+cols[18:]
                    break
            elif int(b[1])-pos>50:
                bicCds='-'
                bicClin='-'
                bicDesign='-'
                newCols=cols[:18]+['','','']+cols[18:]
                break
    # if it is deletion
    elif altAl=='-':
        # we convert it to the BIC format, not increasing value of position
        # because now we use ANNOVAR
        # and make sequence
        delNucs=refAl
        #pos+=len(altAl)
        # make mutated sequence region
        regionStart=pos-200
        regionEnd=pos+200
        altAlleleSeq=refSeqs[brca][regionStart-1:pos-1]+refSeqs[brca][pos+len(delNucs)-1:regionEnd]
        # we search coordinate from BIC that is not less far than 200 bp
        for z,b in enumerate(brca_bic_db[cols[3]][0]):
            # if we found such coordinate
            if abs(pos-int(b[1]))<=50 and b[0]=='deletion':
                # we align our mutated sequence with seq from db
                a=pairwise2.align.globalms(b[8],altAlleleSeq,2,-1,-1.53,0)
                # This alignment has the following format:
                # [('---TCAGTAACAAATGCTCCTATA-', 'GGCTCAGTAACAAATGCTCCTATAA', 21.0, 0, 25),...]
                # we calculate offset for both sides of both sequences
                ls1=len(a[0][0])-len(a[0][0].lstrip('-'))
                ls2=len(a[0][1])-len(a[0][1].lstrip('-'))
                rs1=len(a[0][0].rstrip('-'))
                rs2=len(a[0][1].rstrip('-'))
                if a[0][0][max(ls1,ls2):min(rs1,rs2)]==a[0][1][max(ls1,ls2):min(rs1,rs2)]:
                    bicCheck=True
                    bicCds=b[6]
                    bicClin=b[7]
                    bicDesign=b[8]
                    newCols=cols[:18]+b[6:9]+cols[18:]
                    break
            elif int(b[1])-pos>50:
                bicCds='-'
                bicClin='-'
                bicDesign='-'
                newCols=cols[:18]+['','','']+cols[18:]
                break
    resultFile.write('\t'.join(newCols)+'\n')
    showPercWork(i+1,allWork)
resultFile.close()
print()
