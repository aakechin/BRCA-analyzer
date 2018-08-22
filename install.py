# This script configurates BRCA-analyzer

import argparse,os
import subprocess as sp

thisDir=os.path.dirname(os.path.realpath(__file__))+'/'

par=argparse.ArgumentParser(description='This script configures BRCA-analyzer')
par.add_argument('--bwa','-bwa',dest='bwaDir',type=str,help='destination of BWA. If it can be started with command bwa, type 0 (Default: 0)',default='0',required=True)
par.add_argument('--samtools','-sam',dest='samDir',type=str,help='destination of samtools. If it can be started with command samtools, type 0 (Default: 0)',default='0',required=True)
##par.add_argument('--bcftools','-bcf',dest='bcfDir',type=str,help='destination of bcftools. If it can be started with command bcftools, type 0 (Default: 0)',default='0',required=True)
par.add_argument('--picard','-picard',dest='picardDir',type=str,help='destination of picard.jar (version 2.0.1). For example, ~/picard-2.0.1/dist/',required=True)
par.add_argument('--gatk','-gatk',dest='gatkDir',type=str,help='destination of GenomeAnalysisTK.jar (version 3.6). For example, ~/GenomeAnalysisTK-3.6/',required=True)
par.add_argument('--snpeff','-snpeff',dest='snpeffDir',type=str,help='destination of snpEff.jar. For example, ~/snpEff/',required=True)
par.add_argument('--annovar','-annovar',dest='annovarDir',type=str,help='destination of annovar. For example, ~/annovar/',required=True)
args=par.parse_args()

if args.bwaDir=='0':
    args.bwaDir=''
elif args.bwaDir[-1]!='/':
    args.bwaDir+='/'
if args.samDir=='0':
     args.samDir=''
elif args.samDir[-1]!='/':
    args.samDir+='/'
##if args.bcfDir=='0':
##     args.bcfDir=''
##elif args.bcfDir[-1]!='/':
##    args.bcfDir+='/'
if args.picardDir[-1]!='/':
    args.picardDir+='/'
if args.gatkDir[-1]!='/':
    args.gatkDir+='/'
if args.snpeffDir[-1]!='/':
    args.snpeffDir+='/'
if args.annovarDir[-1]!='/':
    args.annovarDir+='/'

file=open(thisDir+'config.txt','w')
##file.write('\n'.join([args.bwaDir,args.samDir,args.bcfDir,args.picardDir,args.gatkDir,args.snpeffDir,args.annovarDir]))
file.write('\n'.join([args.bwaDir,args.samDir,args.picardDir,args.gatkDir,args.snpeffDir,args.annovarDir]))
file.close()
