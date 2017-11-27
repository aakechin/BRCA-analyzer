# This script recalculates allele distribution through mpileup of samtools

import re,os
import argparse
import subprocess as sp

thisDir=os.path.dirname(os.path.realpath(__file__))+'/'
configs=open(thisDir+'config.txt').read().split('\n')

par=argparse.ArgumentParser(description='This script recalculates allele distribution through mpileup of samtools')
par.add_argument('--input','-in',dest='inputFile',type=str,help='file with avinput',required=True)
par.add_argument('--min-BQ','-bq',dest='minBaseQual',type=int,help='Minimal base quality (BQ) for considering this base (parameter for samtools. Default: 0',default=0)
args=par.parse_args()

file=open(args.inputFile)
# First of all we make BED-file for samtools mpileup
coords=[]
for string in file:
    cols=string[:-1].split('\t')
    coords.append(cols[5:7])
##    bedFile.write('\t'.join(cols[5:7])+'\n')
##bedFile.close()
file.close()
bamFile=args.inputFile[:-len('.unifiedGenotyper.ann.avinput')]+'.bam'
mpileup={}
for coord in coords:
    out=sp.check_output(configs[1]+'samtools mpileup -Q '+str(args.minBaseQual)+' -q 1 -A -t DP -f '+thisDir+'ref/ucsc.hg19.fasta -d 100000 -L 1000000 -r '+':'.join(coord)+'-'+coord[1]+' '+bamFile+' > '+bamFile+'.mpileup0',shell=True,stderr=sp.STDOUT).decode('utf-8')
    # Read firstly mpileup file
    mpileupFile=open(bamFile+'.mpileup0')
    for string2 in mpileupFile:
        cols=string2[:-1].split('\t')
        mpileup['_'.join(cols[0:2])]=cols
inputFile=open(args.inputFile)
resultFile=open(args.inputFile[:-8]+'.ad_recal.avinput','w')
for string in inputFile:
    # Read string of file with avinput
    inputCols=string[:-1].split('\t')
    cols=mpileup['_'.join(inputCols[5:7])]
    if '_'.join(inputCols[5:7]) not in mpileup.keys():
        print('ERROR: Coordinate from inputFile is absent in the mpileup file!')
        print(inputCols[5:7])
        exit(0)
    cols[2]=cols[2].upper()
    alls={'ref':0}
    nucs=['A','T','G','C','a','t','g','c']
    i=0
    p1=re.compile('(?:[ATGCatgc\.\,]\+(\d+)|[ATGCatgc\.\,]\-(\d+)|[ATGCatgc\.\,\*])')
    skip=0
    for m in p1.finditer(cols[4]):
        if skip>0:
            skip-=1
            continue
        start=m.start()
        found=m.group()
        # If we found insertion (+)
        if m.groups()[0]!=None:
            num=int(m.groups()[0])
            al=cols[2]+'>'+found[0].replace('.',cols[2]).replace(',',cols[2]).upper()+cols[4][start+2+len(str(num)):start+2+len(str(num))+num].upper()
            skip=num
        # If we found deletion (-)
        elif m.groups()[1]!=None:
            num=int(m.groups()[1])
            al=cols[2]+cols[4][start+2+len(str(num)):start+2+len(str(num))+num].upper()+'>'+found[0].replace('.',cols[2]).replace(',',cols[2]).upper()
            skip=num
        # If this reference allele
        elif found=='.' or found==',':
            al='ref'
        # If this read has deleted nucleotide in this position
        elif found=='*':
            al=cols[2]+'>-'
        # If this is substitution
        else:
            al=cols[2]+'>'+found.upper()
        if al not in alls.keys():
            alls[al]=1
        else:
            alls[al]+=1
    inputCols14=inputCols[14].split(':')
    oldCov=inputCols[14].split(':')[1].split(',')
    ref=inputCols[8]
    alt=inputCols[9]
    altCov=0
    # When we see as an example CAAAAAAAAAAAAG>CG or CAAAAAAAAAAAAG>CAAAAAAAAAAAG
    if (len(inputCols[9])>1 and inputCols[9][1:] in inputCols[8]) or (len(inputCols[8])>1 and inputCols[8][1:] in inputCols[9]):
        z=1
        while z<min(len(inputCols[8]),len(inputCols[9])):
            if inputCols[8][-z]==inputCols[9][-z]:
                ref=ref[:-1]
                alt=alt[:-1]
            else:
                break
            z+=1
    # When we search C>CTCCATTGCAG in {'ref': 52, 'C>TTCCATTGCAG': 6, 'C>-': 26, 'C>T': 4}
    elif len(inputCols[9])>1 and inputCols[9][0]==inputCols[8][0] and len(inputCols[8])==1:
        for key,item in alls.items():
            if key=='ref': continue
            if '->'+inputCols[9][1:]=='->'+key.split('>')[1][1:]:
                altCov+=item
    # When we search CTCCATTGCAG>C in {'ref': 52, 'CTCCATTGCAG>T': 6, 'C>-': 26, 'C>T': 4}
    elif len(inputCols[8])>1 and inputCols[9][0]==inputCols[8][0] and len(inputCols[9])==1:
        for key,item in alls.items():
            if key=='ref': continue
            if inputCols[8][1:]+'>-'==key.split('>')[0][1:]+'>-':
                altCov+=item
    # When we search A>G or A>AC in {A>GC,A>AC}
    elif len(ref)==1 and len(alt)==1:
        for key,item in alls.items():
            if key=='ref': continue
            if ref+'>'+alt==key.split('>')[0][0]+'>'+key.split('>')[1][0]:
                altCov+=item
    if ref+'>'+alt not in alls.keys() and altCov==0:
        print('ERROR: There is no such mutation in BAM-file!')
        print(args.inputFile)
        print('Coordinate:','_'.join(inputCols[5:7]))
        print('Mutation from variation file:',inputCols[8]+'>'+inputCols[9])
        print('Mutations from BAM-file:',alls)
        exit(0)
    elif altCov>0:
        resultFile.write('\t'.join(inputCols[:14])+'\t'+inputCols14[0]+':'+str(alls['ref'])+','+str(altCov)+':'+':'.join(inputCols14[2:])+'\n')
    else:
        resultFile.write('\t'.join(inputCols[:14])+'\t'+inputCols14[0]+':'+str(alls['ref'])+','+str(alls[ref+'>'+alt])+':'+':'.join(inputCols14[2:])+'\n')
##    print(inputCols[5:8],inputCols[8:10],oldCov,alls)
inputFile.close()
mpileupFile.close()
resultFile.close()
    
