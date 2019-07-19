# This script do all processing of BRCA data

# Section of importing modules
import argparse
import subprocess as sp
import os
import sys
from multiprocessing import Pool,Queue
import multiprocessing as mp
import glob
import re,math

# Global variables
global thisDir,configs,allWork,threads,args
thisDir=os.path.dirname(os.path.realpath(__file__))+'/'
configs=open(thisDir+'config.txt').read().split('\n')

# Section of functions
def showPercWork(done,allWork):
    percDoneWork=round((done/allWork)*100,2)
    sys.stdout.write("\r"+str(percDoneWork)+"%")
    sys.stdout.flush()

def processPatient(inputArgs):
    varCaller='pisces'
    readsFile1=inputArgs[0]
    readsFile2=inputArgs[1]
    outDir=inputArgs[2]
    patNum=inputArgs[3]
    if not os.path.isdir(outDir+'patient_'+patNum):
        os.mkdir(outDir+'patient_'+patNum)
    # Map reads to the reference
    if (not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+'.sam')
        and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+'.sam.gz')
        and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+'.bam')
        and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+'.bam.gz')):
        if readsFile2:
            output=sp.check_output(configs[0]+"bwa mem "+configs[7]+''
                                   ' '+readsFile1+' '+readsFile2+' -t '+threads+' | gzip > '+outDir+''
                                   'patient_'+patNum+'/patient_'+patNum+'.sam.gz',
                                   shell=True,stderr=sp.STDOUT)
        else:
            output=sp.check_output(configs[0]+"bwa mem "+configs[7]+''
                                   ' '+readsFile1+' -t '+threads+' | gzip > '+outDir+''
                                   'patient_'+patNum+'/patient_'+patNum+'.sam.gz',
                                   shell=True,stderr=sp.STDOUT)
        if 'fail to locate the index files' in str(output):
            print('#'*10,'\nERROR: File with reference sequence was not found!\n'
                  'It should be:',thisDir+'ref/ucsc.hg19.fasta\n','#'*10)
            exit(1)
    # Add read groups
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.bam') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.bam.gz'):
        output=sp.check_output("java -jar "+configs[2]+"picard.jar "
                               "AddOrReplaceReadGroups INPUT="+outDir+'patient_'+patNum+'/patient_'+patNum+'.sam.gz'
                               ' OUTPUT='+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.bam'
                               ' SORT_ORDER=coordinate CREATE_INDEX=TRUE RGLB=MiSeq RGPL=Illumina'
                               ' RGPU=barcode RGSM=patient_'+patNum,
                               shell=True,stderr=sp.STDOUT)
    # Create targets for realignment
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.intervals') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.intervals.gz'):
        output=sp.check_output("java -jar "+configs[3]+"GenomeAnalysisTK.jar"
                               ' -T RealignerTargetCreator'
                               ' -R '+configs[7]+''
                               ' -I '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.bam'
                               ' -o '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.intervals'
                               ' -L chr13:32889617-32973809'
                               ' -L chr17:41196312-41279700 -nt '+threads,
                               shell=True,stderr=sp.STDOUT)
    # Realign
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.bam') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.bam.gz'):
        output=sp.check_output("java -jar "+configs[3]+"GenomeAnalysisTK.jar"
                               ' -T IndelRealigner'
                               ' -R '+configs[7]+''
                               ' -I '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.bam'
                               ' -o '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.bam'
                               ' -targetIntervals '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.intervals'
                               ' -L chr13:32889617-32973809'
                               ' -L chr17:41196312-41279700',
                               shell=True,stderr=sp.STDOUT)
    # Recalibration
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal_data.table') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal_data.table.gz'):
        output=sp.check_output("java -jar "+configs[3]+"/GenomeAnalysisTK.jar"
                               ' -T BaseRecalibrator'
                               ' -R '+configs[7]+''
                               ' -I '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.bam'
                               ' -o '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal_data.table'
                               ' -knownSites '+thisDir+'annotation_databases/clinvar_20160831.vcf.gz'
                               ' -L chr13:32889617-32973809'
                               ' -L chr17:41196312-41279700 -nct '+threads,
                               shell=True,stderr=sp.STDOUT)
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.post_recal_data.table') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.post_recal_data.table.gz'):
        output=sp.check_output("java -jar "+configs[3]+"GenomeAnalysisTK.jar"
                               ' -T BaseRecalibrator'
                               ' -R '+configs[7]+''
                               ' -I '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.bam'
                               ' -o '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.post_recal_data.table'
                               ' -knownSites '+thisDir+'annotation_databases/clinvar_20160831.vcf.gz'
                               ' -BQSR '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal_data.table'
                               ' -L chr13:32889617-32973809'
                               ' -L chr17:41196312-41279700 -nct '+threads,
                               shell=True,stderr=sp.STDOUT)
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal.bam') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal.bam.gz'):
        output=sp.check_output("java -jar "+configs[3]+"GenomeAnalysisTK.jar"
                               ' -T PrintReads'
                               ' -R '+configs[7]+''
                               ' -I '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.bam'
                               ' -o '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.bam'
                               ' -BQSR '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal_data.table'
                               ' -L chr13:32889617-32973809'
                               ' -L chr17:41196312-41279700 -nct '+threads,
                               shell=True,stderr=sp.STDOUT)
    # Cut primer sequences from reads
    cutPrimers=''
    if args.primersFileR1_5 and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                                                   '.sorted.read_groups.realigned.recal.trimmed.sorted.bam'):
        cmd='python3 '+configs[6]+'cutPrimers.py -bam '
        cmd+=outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.bam'
        cmd+=' -pr15 '+args.primersFileR1_5
        if args.primersFileR1_3:
            cmd+=' -pr13 '+args.primersFileR1_3
        if args.primersFileR2_5:
            cmd+=' -pr25 '+args.primersFileR2_5
        if args.primersFileR2_3:
            cmd+=' -pr23 '+args.primersFileR2_3
        cmd+=' -outbam '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.trimmed.bam'
        cmd+=' -outbam2 '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.untrimmed.bam'
        cmd+=' -coord '+args.coordsFile
        cmd+=' -t '+str(args.toolThreads)
        cmd+=' -e '+str(args.errNumber)
        cmd+=' -plb '+str(args.primerLocBuf)
        cmd+=' -minlen '+str(args.minReadLen)
        if args.primer3absent:
            cmd+=' -primer3'
        output=sp.check_output(cmd,shell=True,stderr=sp.STDOUT)
        cutPrimers='trimmed.sorted.'
    elif args.primersFileR1_5 and os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                                                 '.sorted.read_groups.realigned.recal.trimmed.sorted.bam'):
        cutPrimers='trimmed.sorted.'
    # Call variants
    if varCaller=='unifiedGenotyper':
        if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                              '.sorted.read_groups.realigned.recal.'+cutPrimers+'unifiedGenotyper.vcf') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                              '.sorted.read_groups.realigned.recal.'+cutPrimers+'unifiedGenotyper.vcf.gz'):
            cmd=["java -jar "+configs[3]+"GenomeAnalysisTK.jar",
                 ' -T UnifiedGenotyper',' -R '+configs[7],
                 ' -I '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.'+cutPrimers+'bam',
                 ' -o '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.'+cutPrimers+'unifiedGenotyper.vcf',
                 ' -dfrac 1',' -glm BOTH',' -minIndelFrac 0.01',' -mbq 10',' -L chr13:32889617-32973809',
                 ' -L chr17:41196312-41279700 -nt '+threads]
            output=sp.check_output(' '.join(cmd),shell=True,stderr=sp.STDOUT)
    elif varCaller=='freebayes':
        if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                              '.sorted.read_groups.realigned.recal.'+cutPrimers+'freebayes.vcf') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                              '.sorted.read_groups.realigned.recal.'+cutPrimers+'freebayes.vcf.gz'):
            cmd=["/home/andrey/Downloads/freebayes/bin/freebayes",
                 ' -f '+configs[7],
                 '-p 4','-0','-F 0.01','-k',
                 '-r chr13:32889617-32973809','-r chr17:41196312-41279700']
            cmd.append(outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.'+cutPrimers+'bam')
            cmd.append('| /home/andrey/Downloads/freebayes/vcflib/bin/vcfallelicprimitives -kg')
            cmd.append('> '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.'+cutPrimers+'freebayes.vcf')
            output=sp.check_output(' '.join(cmd),shell=True,stderr=sp.STDOUT)
    elif varCaller=='pisces':
        varCaller='genome'
        if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                              '.sorted.read_groups.realigned.recal.'+cutPrimers+varCaller+'.vcf') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                              '.sorted.read_groups.realigned.recal.'+cutPrimers+varCaller+'.vcf.gz'):
            cmd=["dotnet /home/andrey/Downloads/Pisces-master/binaries/5.2.9.122/Pisces_5.2.9.122/Pisces.dll",
##                 '--coveragemethod "exact"','--sbfilter 0.1',
                 '--gvcf false','--minbq 10','--minmq 10','--minvf 0.01',
                 '--minvq 1',
                 '-bam',outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.'+cutPrimers+'bam',
                 '-g /srv/reference_genomes/human/']
            if args.coordsFile:
                cmd.append('-i '+args.coordsFile)
            output=sp.check_output(' '.join(cmd),shell=True,stderr=sp.STDOUT)
            os.rename(outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.'+cutPrimers+'vcf',
                      outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.'+cutPrimers+varCaller+'.vcf')
   # Annotate variants
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal.'+cutPrimers+varCaller+'.ann.vcf') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal.'+cutPrimers+varCaller+'.ann.vcf.gz'):
        output=sp.check_output("java -jar "+configs[4]+"snpEff.jar hg19"
                               ' '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.'+cutPrimers+varCaller+'.vcf'
                               ' > '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.'+cutPrimers+varCaller+'.ann.vcf',
                               shell=True,stderr=sp.STDOUT)
    # Convert VCF to ANNOVAR input file
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal.'+cutPrimers+varCaller+'.ann.avinput') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal.'+cutPrimers+varCaller+'.ann.avinput.gz'):
        cmd=[configs[5]+"convert2annovar.pl",
             ' -format vcf4 -includeinfo','-allsample','-withfreq',
             ' '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.'+cutPrimers+varCaller+'.ann.vcf',
             ' >'+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.'+cutPrimers+varCaller+'.ann.avinput']
        output=sp.check_output(' '.join(cmd),shell=True,stderr=sp.STDOUT)
##    # Recalculate allele distribution
##    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
##                          '.sorted.read_groups.'+cutPrimers+'realigned.recal.unifiedGenotyper.ann.ad_recal.avinput') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
##                          '.sorted.read_groups.'+cutPrimers+'realigned.recal.unifiedGenotyper.ann.ad_recal.avinput.gz'):
##        output=sp.check_output('python3 '+thisDir+'recalculateAlleleDistribution.py'
##                               ' -in '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.'+cutPrimers+'realigned.recal.unifiedGenotyper.ann.avinput',
##                               shell=True,stderr=sp.STDOUT).decode('utf-8')
##        if 'ERROR' in output:
##            print(output)
##            exit(1)
    q1.put(1)
    showPercWork(q1.qsize(),allWork)

# Section of input arguments
par=argparse.ArgumentParser(description='This script do all processing of BRCA sequencing data')
par.add_argument('--readsFiles_r1','-r1',dest='readsFiles1',type=str,help='regular expression for choosing files with R1 reads',required=True)
par.add_argument('--readsFiles_r2','-r2',dest='readsFiles2',type=str,help='regular expression for choosing files with R2 reads (optional)',required=False)
par.add_argument('--readsFiles_N','-rN',dest='readsFilesN',type=str,help='regular expression for choosing files with native R1 reads (including Undetermined) to evaluate effectivity of trimming reads  (alternative). If you do not have trimmed reads, do not use this parameter',required=False)
par.add_argument('--primersFileR1_5','-pr15',dest='primersFileR1_5',type=str,
                 help='fasta-file with sequences of primers on the 5\'-end of R1 reads. '
                 'Use it, only if you need to cut primer sequences from reads',required=False)
par.add_argument('--primersFileR2_5','-pr25',dest='primersFileR2_5',type=str,
                 help='fasta-file with sequences of primers on the 5\'-end of R2 reads. '
                 'Use it, only if you need to cut primer sequences from reads. '
                 'Also, do not use this parameter if you have single-end reads',required=False)
par.add_argument('--primersFileR1_3','-pr13',dest='primersFileR1_3',type=str,
                 help='fasta-file with sequences of primers on the 3\'-end of R1 reads. '
                 'Use it, only if you need to cut primer sequences from reads. '
                 'It is not required. But if it is determined, -pr23 is necessary',required=False)
par.add_argument('--primersFileR2_3','-pr23',dest='primersFileR2_3',type=str,
                 help='fasta-file with sequences of primers on the 3\'-end of R2 reads. '
                 'Use it, only if you need to cut primer sequences from reads. ',required=False)
par.add_argument('--coordinates-file','-coord',dest='coordsFile',type=str,
                 help='file with coordinates of amplicons in the BED-format '
                 '(without column names and locations of primers): chromosome | start | end. '
                 'It is necessary for cutting primer sequences from BAM-file. '
                 'Its order should be the same as for files with primer sequences',required=False)
par.add_argument('--error-number','-err',dest='errNumber',type=int,
                 help='number of errors (substitutions, insertions, deletions) '
                 'that allowed during searching primer sequence in a read sequence. Default: 3',default=3)
par.add_argument('--primer-location-buffer','-plb',dest='primerLocBuf',type=int,
                 help='Buffer of primer location in the read from the start or end of read. '
                 'If this value is zero, than cutPrimers will search for primer sequence '
                 'in the region of the longest primer length. Default: 10',default=10)
par.add_argument('--minimal-read-length','-minlen',dest='minReadLen',type=int,
                 help='minimal length of read after trimming. Default: 30',default=30)
par.add_argument('--primer3-absent','-primer3',dest='primer3absent',action='store_true',
                 help="if primer at the 3'-end may be absent, use this parameter")
par.add_argument('--patientsTable','-pat',dest='patientsTable',type=str,help='table with information about each patient: ngs_num patient_id barcode1 barcode2',required=False)
par.add_argument('--primersCoords','-primer',dest='primersCoords',type=str,help='table with information about amplicon coordinates without column headers: amplicon_number | chromosome | start | end. (Is not required)',required=False)
par.add_argument('--outDir','-out',dest='outDir',type=str,help='directory for output',required=True)
par.add_argument('--threads','-th',dest='threads',type=int,help='number of threads',default=2)
par.add_argument('--tool-threads','-tt',dest='toolThreads',type=int,help='number of threads for each tool. Number of --threads multiplied by the number of --tool-threads must not exceed number of CPU cores',default=1)
par.add_argument('--run-name','-run',dest='runName',type=str,help='Name of run. This name will be added into the name of the final output file. Default: BRCA',required=False,default='BRCA')
par.add_argument('--without-joinment','-notjoin',dest='notToJoin',action='store_true',help='use this parameter if you only want to process patients reads separately without joining them (useful if you have an opportunity to separate processing onto several machines)')
par.add_argument('--only-join','-onlyjoin',dest='onlyJoin',action='store_true',help='use this parameter if you only want to join already processed patients reads')
par.add_argument('--language','-lang',dest='lang',type=str,help='Language of report and text on figures (russian or english). Default: english',default='english')
args=par.parse_args()

varCaller='pisces'
if varCaller=='pisces':
    varCaller='genome'
langs=['russian','english']
if args.lang=='ru': args.lang='russian'
elif args.lang=='en': args.lang='english'
if args.lang not in langs:
    print('#'*10,'\nWARNING! Chosen language is not accepted. Use default english...')
if args.notToJoin and args.onlyJoin:
    print("ERROR: only one of parameters can be used, -notjoin or -onlyjoin. But you have tried to use both")
    exit(1)
readsFiles1=sorted(glob.glob(args.readsFiles1))
if len(readsFiles1)==0:
    print('ERROR: no files for reads R1 were selected!'); exit(1)
if args.readsFiles2:
    readsFiles2=sorted(glob.glob(args.readsFiles2))
    if len(readsFiles2)==0:
        print('ERROR: no files for reads R2 were selected!'); exit(1)
else:
    readsFiles2=[False]*len(readsFiles1)
if args.readsFilesN:
    readsFilesN=sorted(glob.glob(args.readsFilesN))
    if len(readsFilesN)==0:
        print('ERROR: no files for native reads were selected!'); exit(1)
if args.patientsTable:
    patientsTable=args.patientsTable
    if not os.path.exists(patientsTable):
        print('ERROR: patients table file does not exist!'); exit(1)
    else:
        patientsTable=os.path.abspath(patientsTable)
if args.primersCoords and not os.path.exists(args.primersCoords):
    print('ERROR: primers coords table file does not exist!'); exit(1)
elif not args.primersCoords:
    args.primersCoords=thisDir+'primers_coords_default.csv'
outDir=args.outDir
if outDir[-1]!='/': outDir+='/'
# Create output directory
if not os.path.isdir(outDir):
    os.mkdir(outDir)
t=int(args.threads)
if not args.onlyJoin:
    rect=int(mp.cpu_count()/args.toolThreads)
    if t>rect:
        print('WARNING: number of threads ('+str(t)+') is too high! Recomended number of threads is '+str(rect))
        print('Continue anyway? (y/n)')
        text=input()
        if text!='Y' and text!='y':
            exit(1)
    threads=str(args.toolThreads)
    output=sp.check_output('df '+outDir,shell=True,stderr=sp.STDOUT).decode('utf-8')
    available=int(output.split('\n')[1].split()[3])
    neededFreeSpace=0
    for r1,r2 in zip(readsFiles1,readsFiles2):
        neededFreeSpace+=os.path.getsize(r1)+os.path.getsize(r2)
    neededFreeSpace=int(neededFreeSpace*4/1024)
    if available<neededFreeSpace:
        print('WARNING: there is not enough free space for analysis. You need about '+str(neededFreeSpace)+' Kb. But available free space is '+str(available)+' Kb')
        print('Continue anyway? (y/n)')
        text=input()
        if text!='Y' and text!='y':
            exit(1)
    # Check that file with reference genome sequence exists
    if not os.path.exists(configs[7]):
        print('#'*10,'\nERROR: File with reference sequence was not found!\n'
                  'It should be:',thisDir+'ref/ucsc.hg19.fasta\n'+'#'*10)
        exit(1)

# Check that file with reference sequence for chromosome 13 and 17 exists
if not os.path.exists(thisDir+'ref/human_g1k_v37_chr13+17.fasta'):
    if os.path.exists(thisDir+'ref/human_g1k_v37_chr13+17.fasta.tar.gz'):
        print('#'*10,'\nERROR: You forgot to untar archive with sequences of chromosomes 13 and 17 in the "ref" directory!\n'+'#'*10)
    else:
        print('#'*10,'\nERROR: Sequences of chromosomes 13 and 17 is absent in the "ref" directory!\n'
              'It should look like:',thisDir+'ref/human_g1k_v37_chr13+17.fasta\n'+'#'*10)
    exit(1)

# Section of procedures
# reads quality evaluation and cutting primers are distinct processes
# that are done separately
# Extract numbers of patients
# For doing it we compare all numbers from readsFiles names and
# save different ones
if len(readsFiles1)>1:
    p=re.compile('\d+')
    indexP=re.compile('_\d+_\d+')
    readsFiles1Name1=readsFiles1[0][readsFiles1[0].rfind('/')+1:]
    readsFiles1Name2=readsFiles1[1][readsFiles1[1].rfind('/')+1:]
    m1=list(set(p.findall(readsFiles1Name1)))
    m2=list(set(p.findall(readsFiles1Name2)))
    nums=set(m1).intersection(m2)
    patNums=[]
    for rf in readsFiles1:
        rf=rf[rf.rfind('/')+1:]
        m=p.findall(rf)
        m2=indexP.findall(rf)
        if len(m)>0 and len(m2)==0:
            for n in nums:
                try:
                    m.remove(n)
                except ValueError:
                    print('INTERNAL ERROR: There is no such element "'+n+'" in nums:')
                    print(readsFiles1[0])
                    print(nums)
                    exit(1)
            patNums.append(m[0])
        # We have index numbers in our reads files names
        elif len(m)>0 and len(m2)>0:
            patNums.append(m2[0])
        else:
            print('INTERNAL ERROR: There is no numbers in the names of reads files!')
            exit(1)
else:
    patNums=['1']
allWork=len(patNums)
outDirs=[outDir]*allWork
if not args.onlyJoin:
    # So on this step we start from mapping reads
    q1=Queue()
    p=Pool(t)
    print('Process reads...')
    showPercWork(0,allWork)
    p.map(processPatient,zip(readsFiles1,readsFiles2,outDirs,patNums))
    print('\nDone')
if args.notToJoin:
    exit(1)
elif args.onlyJoin:
    ds=glob.glob(outDir+'patient_*/*.unifiedGenotyper.ann.avinput')
    if len(ds)<len(patNums):
        print('#'*10)
        print('WARNING: Number of files with completely processed reads is less than number of patients!')
        print('Number of patients:',len(patNums))
        print('Number of files with completely processed reads (*.unifiedGenotyper.ann.avinput):',len(ds))
        print('#'*10)
# Join mutations of all patients to one file
print('Joining mutations of all patients...')
output=sp.check_output("python3 "+thisDir+"joinMutations.py -v '"+outDir+"patient_*/*."+varCaller+".ann.avinput' -o "+outDir+'allPatients.'+varCaller+'.ann.avinput',shell=True,stderr=sp.STDOUT)
if 'ERROR' in str(output):
    print(output.decode('utf-8'))
    exit(1)
print('Done')
# Annotate them with ANNOVAR
print('Annotating mutations of all patients with ANNOVAR...')
output=sp.check_output(configs[5]+"table_annovar.pl "+outDir+'allPatients.'+varCaller+'.ann.avinput '+configs[5]+'humandb/ -buildver hg19 -out '+outDir+'allPatients.'+varCaller+'.ann.avinput.ann -remove -protocol refGene,cosmic70,esp6500siv2_all,exac03,kaviar_20150923,1000g2015aug_all,avsnp150,clinvar_20180603,ljb26_all -operation g,f,f,f,f,f,f,f,f -nastring . -otherinfo',shell=True,stderr=sp.STDOUT)
if 'ERROR' in str(output):
    print(output.decode('utf-8'))
    exit(1)
print('Done')
# Add positions in old references
print('Adding positions in old references...')
if args.patientsTable:
    output=sp.check_output("python3 "+thisDir+"addPosOldRef.py -in "+outDir+'allPatients.'+varCaller+'.ann.avinput.ann.hg19_multianno.txt -pt '+patientsTable+' -pl '+'_'.join(patNums),shell=True,stderr=sp.STDOUT)
else:
    output=sp.check_output("python3 "+thisDir+"addPosOldRef.py -in "+outDir+'allPatients.'+varCaller+'.ann.avinput.ann.hg19_multianno.txt -pl '+'_'.join(patNums),shell=True,stderr=sp.STDOUT)
if 'ERROR' in str(output):
    print(output.decode('utf-8'))
    exit(1)
print('Done')
# Annotate by Clinical information
print('Annotating with BIC...')
output=sp.check_output("python3 "+thisDir+"annotateClinVars.py -i "+outDir+'allPatients.'+varCaller+'.ann.avinput.ann.hg19_multianno.withOurCoordinates.xls',shell=True,stderr=sp.STDOUT)
if 'ERROR' in str(output):
    print(output.decode('utf-8'))
    exit(1)
print('Done')
# Convert to Excel-format
print('Converting to Excel-format...')
output=sp.check_output("python3 "+thisDir+"convertResultToExcel.py"
                       ' -i '+outDir+'allPatients.'+varCaller+'.ann.avinput.ann.hg19_multianno.withOurCoordinates.clinical.xls'
                       ' -o '+outDir+'allPatients.'+args.runName+'.'+varCaller+'.ann.avinput.ann.hg19_multianno.withOurCoordinates.clinical.excel.xls',
                       shell=True,stderr=sp.STDOUT)
if 'ERROR' in str(output):
    print(output.decode('utf-8'))
    exit(1)
print('Done')
if args.primersFileR1_5:
    cutPrimers='trimmed.sorted.'
else:
    cutPrimers=''
if args.patientsTable:
    # Evaluating coverage
    print('Evaluating coverage...')
    process=sp.Popen(["python3",thisDir+"countCovForAmplic.py","-mmm","min",
                      '-coord',args.primersCoords,'-pat',patientsTable,'-bam',
                      ''+outDir+'patient_*/*.sorted.read_groups.realigned.recal.'+cutPrimers+'bam','-res',
                      outDir+args.runName+'_amplicon_coverage_min.xls','-th',str(t*int(args.toolThreads)),'-pl','_'.join(patNums)],
                     stdout=sp.PIPE,shell=False,universal_newlines=True)
    for c in iter(process.stdout.readline,''):
        sys.stdout.write("\r"+str(c).replace('\n',''))
        sys.stdout.flush()
    print()
    process=sp.Popen(["python3",thisDir+"countCovForAmplic.py","-mmm","mean",
                      '-coord',args.primersCoords,'-pat',patientsTable,'-bam',
                      ''+outDir+'patient_*/*.sorted.read_groups.realigned.recal.'+cutPrimers+'bam','-res',
                      outDir+args.runName+'_amplicon_coverage_mean.xls','-th',str(t*int(args.toolThreads)),'-pl','_'.join(patNums)],
                     stdout=sp.PIPE,shell=False,universal_newlines=True)
    for c in iter(process.stdout.readline,''):
        sys.stdout.write("\r"+str(c).replace('\n',''))
        sys.stdout.flush()
    print()
else:
    # Evaluating coverage
    print('Evaluating coverage...')
    process=sp.Popen(["python3",thisDir+"countCovForAmplic.py","-mmm","min",
                      '-coord',args.primersCoords,'-bam',
                      ''+outDir+'patient_*/*.sorted.read_groups.realigned.recal.'+cutPrimers+'bam','-res',
                      outDir+args.runName+'_amplicon_coverage_min.xls','-th',str(t*int(args.toolThreads)),'-pl','_'.join(patNums)],
                     stdout=sp.PIPE,shell=False,universal_newlines=True)
    for c in iter(process.stdout.readline,''):
        sys.stdout.write("\r"+str(c).replace('\n',''))
        sys.stdout.flush()
    print()
    process=sp.Popen(["python3",thisDir+"countCovForAmplic.py","-mmm","mean",
                      '-coord',args.primersCoords,'-bam',
                      ''+outDir+'patient_*/*.sorted.read_groups.realigned.recal.'+cutPrimers+'bam','-res',
                      outDir+args.runName+'_amplicon_coverage_mean.xls','-th',str(t*int(args.toolThreads)),'-pl','_'.join(patNums)],
                     stdout=sp.PIPE,shell=False,universal_newlines=True)
    for c in iter(process.stdout.readline,''):
        sys.stdout.write("\r"+str(c).replace('\n',''))
        sys.stdout.flush()
    print()
# Evaluate read statistics
print('Evaluating reads...')
if args.readsFilesN:
    if args.patientsTable:
        output=sp.check_output('python3 '+thisDir+'getPercentOfProperlyTrimmedReads.py '
                               '-nat "'+args.readsFilesN+'" -trim "'+args.readsFiles1+'" -pat '+patientsTable+' -out '+outDir+args.runName+'.reads_statistics.xls',shell=True,stderr=sp.STDOUT)
    else:
        output=sp.check_output('python3 '+thisDir+'getPercentOfProperlyTrimmedReads.py '
                               '-nat "'+args.readsFilesN+'" -trim "'+args.readsFiles1+'" -out '+outDir+args.runName+'.reads_statistics.xls',shell=True,stderr=sp.STDOUT)
elif args.primersFileR1_5:
    if args.patientsTable:
        output=sp.check_output(' '.join(['python3',
                                         thisDir+'getPercentOfProperlyTrimmedReads.py',
                                         '-nat','"'+args.readsFiles1+'"',
                                         '-trim',
                                         '"'+outDir+'patient_*/*.sorted.read_groups.realigned.recal.trimmed.sorted.bam'+'"',
                                         '-pat',args.patientsTable,
                                         '-out',outDir+args.runName+'.reads_statistics.xls']),
                               shell=True,stderr=sp.STDOUT)
    else:
        output=sp.check_output(' '.join(['python3',
                                         thisDir+'getPercentOfProperlyTrimmedReads.py',
                                         '-nat','"'+args.readsFiles1+'"',
                                         '-trim',
                                         '"'+outDir+'patient_*/*.sorted.read_groups.realigned.recal.trimmed.sorted.bam'+'"',
                                         '-out',outDir+args.runName+'.reads_statistics.xls']),
                               shell=True,stderr=sp.STDOUT)
else:
    if args.patientsTable:
        output=sp.check_output('python3 '+thisDir+'getPercentOfProperlyTrimmedReads.py '
                           '-nat "'+args.readsFiles1+'" -pat '+patientsTable+' -out '+outDir+args.runName+'.reads_statistics.xls',shell=True,stderr=sp.STDOUT)
    else:
        output=sp.check_output('python3 '+thisDir+'getPercentOfProperlyTrimmedReads.py '
                           '-nat "'+args.readsFiles1+'" -out '+outDir+args.runName+'.reads_statistics.xls',shell=True,stderr=sp.STDOUT)
print('Done')
if args.primersCoords:
    print('Evaluating uniformity of coverage...')
    output=sp.check_output('python3 '+thisDir+'drawUniformityFigure.py '
                           '-cov '+outDir+args.runName+'_amplicon_coverage_mean.xls '
                           '-out '+outDir+args.runName+'_amplicon_coverage_mean.uniformity.tiff '
                           '-lang '+args.lang,shell=True,stderr=sp.STDOUT)
print('Done')
print('Creating report...')
if args.primersCoords:
    output=sp.check_output('python3 '+thisDir+'makeReport.py '
                            '-read '+outDir+args.runName+'.reads_statistics.xls '
                           '-cov '+outDir+args.runName+'_amplicon_coverage_mean.xls '
                           '-res '+outDir+'allPatients.'+args.runName+'.'+varCaller+'.ann.avinput.ann.hg19_multianno.withOurCoordinates.clinical.excel.xls '
                           '-fig '+outDir+args.runName+'_amplicon_coverage_mean.uniformity.tiff '
                           '-out '+outDir+args.runName+'.report.docx '
                           '-lang '+args.lang,shell=True,stderr=sp.STDOUT)
else:
    output=sp.check_output('python3 '+thisDir+'makeReport.py '
                            '-read '+outDir+args.runName+'.reads_statistics.xls '
                           '-res '+outDir+'allPatients.'+args.runName+'.'+varCaller+'.ann.avinput.ann.hg19_multianno.withOurCoordinates.clinical.excel.xls '
                           '-out '+outDir+args.runName+'.report.docx '
                           '-lang '+args.lang,shell=True,stderr=sp.STDOUT)
print('Done')

