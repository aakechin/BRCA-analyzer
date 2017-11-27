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
global thisDir,configs,allWork,threads
thisDir=os.path.dirname(os.path.realpath(__file__))+'/'
configs=open(thisDir+'config.txt').read().split('\n')

# Section of functions
def showPercWork(done,allWork):
    percDoneWork=round((done/allWork)*100,2)
    sys.stdout.write("\r"+str(percDoneWork)+"%")
    sys.stdout.flush()

def processPatient(inputArgs):
    readsFile1=inputArgs[0]
    readsFile2=inputArgs[1]
    outDir=inputArgs[2]
    patNum=inputArgs[3]
    if inputArgs[4]:
        fixMis=' -fixMisencodedQuals'
    else:
        fixMis=''
    if not os.path.isdir(outDir+'patient_'+patNum):
        os.mkdir(outDir+'patient_'+patNum)
    # Map reads to the reference
    if (not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+'.sam')
        and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+'.sam.gz')
        and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+'.bam')
        and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+'.bam.gz')):
        if readsFile2:
            output=sp.check_output(configs[0]+"bwa mem "+thisDir+'ref/ucsc.hg19.fasta'
                                   ' '+readsFile1+' '+readsFile2+' -t '+threads+' | gzip > '+outDir+''
                                   'patient_'+patNum+'/patient_'+patNum+'.sam.gz',
                                   shell=True,stderr=sp.STDOUT)
        else:
            output=sp.check_output(configs[0]+"bwa mem "+thisDir+'ref/ucsc.hg19.fasta'
                                   ' '+readsFile1+' -t '+threads+' | gzip > '+outDir+''
                                   'patient_'+patNum+'/patient_'+patNum+'.sam.gz',
                                   shell=True,stderr=sp.STDOUT)
    # Convert SAM-file to BAM
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+'.bam') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+'.bam.gz'):
        output=sp.check_output(configs[1]+"samtools view -b "+outDir+'patient_'+patNum+'/patient_'+patNum+'.sam.gz '
                               '-o '+outDir+'patient_'+patNum+'/patient_'+patNum+'.bam'
                               ' -@ '+str(int(threads)-1),
                               shell=True,stderr=sp.STDOUT)
    # Sort BAM-file
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.bam') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.bam.gz'):
        output=sp.check_output(configs[1]+"samtools sort -T /tmp/aln_"+patNum+".sorted "
                               "-o "+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.bam'
                               ' '+outDir+'patient_'+patNum+'/patient_'+patNum+'.bam'
                               ' -@ '+str(int(threads)-1),
                               shell=True,stderr=sp.STDOUT)
    # Add read groups
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.bam') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.bam.gz'):
        output=sp.check_output("java -jar "+configs[2]+"picard.jar "
                               "AddOrReplaceReadGroups INPUT="+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.bam'
                               ' OUTPUT='+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.bam'
                               ' SORT_ORDER=coordinate RGLB=MiSeq RGPL=Illumina RGPU=barcode RGSM=patient_'+patNum,
                               shell=True,stderr=sp.STDOUT)
    # Index BAM-file
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.bam.bai') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.bam.bai.gz'):
        output=sp.check_output(configs[1]+"samtools index "+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.bam',
                               shell=True,stderr=sp.STDOUT)
    # Create targets for realignment
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.intervals') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.intervals.gz'):
        output=sp.check_output("java -jar "+configs[3]+"GenomeAnalysisTK.jar"
                               ' -T RealignerTargetCreator'
                               ' -R '+thisDir+'ref/ucsc.hg19.fasta'
                               ' -I '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.bam'
                               ' -o '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.intervals'
                               ' -L chr13:32889617-32973809'
                               ' -L chr17:41196312-41277500'+fixMis+' -nt '+threads,
                               shell=True,stderr=sp.STDOUT)
    # Realign
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.bam') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.bam.gz'):
        output=sp.check_output("java -jar "+configs[3]+"GenomeAnalysisTK.jar"
                               ' -T IndelRealigner'
                               ' -R '+thisDir+'ref/ucsc.hg19.fasta'
                               ' -I '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.bam'
                               ' -o '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.bam'
                               ' -targetIntervals '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.intervals'
                               ' -L chr13:32889617-32973809'
                               ' -L chr17:41196312-41277500'+fixMis,
                               shell=True,stderr=sp.STDOUT)
    # Recalibration
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal_data.table') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal_data.table.gz'):
        output=sp.check_output("java -jar "+configs[3]+"/GenomeAnalysisTK.jar"
                               ' -T BaseRecalibrator'
                               ' -R '+thisDir+'ref/ucsc.hg19.fasta'
                               ' -I '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.bam'
                               ' -o '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal_data.table'
                               ' -knownSites '+thisDir+'annotation_databases/clinvar_20160831.vcf.gz'
                               ' -L chr13:32889617-32973809'
                               ' -L chr17:41196312-41277500 -nct '+threads,
                               shell=True,stderr=sp.STDOUT)
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.post_recal_data.table') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.post_recal_data.table.gz'):
        output=sp.check_output("java -jar "+configs[3]+"GenomeAnalysisTK.jar"
                               ' -T BaseRecalibrator'
                               ' -R '+thisDir+'ref/ucsc.hg19.fasta'
                               ' -I '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.bam'
                               ' -o '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.post_recal_data.table'
                               ' -knownSites '+thisDir+'annotation_databases/clinvar_20160831.vcf.gz'
                               ' -BQSR '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal_data.table'
                               ' -L chr13:32889617-32973809'
                               ' -L chr17:41196312-41277500 -nct '+threads,
                               shell=True,stderr=sp.STDOUT)
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal.bam') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal.bam.gz'):
        output=sp.check_output("java -jar "+configs[3]+"GenomeAnalysisTK.jar"
                               ' -T PrintReads'
                               ' -R '+thisDir+'ref/ucsc.hg19.fasta'
                               ' -I '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.bam'
                               ' -o '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.bam'
                               ' -BQSR '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal_data.table'
                               ' -L chr13:32889617-32973809'
                               ' -L chr17:41196312-41277500 -nct '+threads,
                               shell=True,stderr=sp.STDOUT)
    # Call variants
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal.unifiedGenotyper.vcf') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal.unifiedGenotyper.vcf.gz'):
        output=sp.check_output("java -jar "+configs[3]+"GenomeAnalysisTK.jar"
                               ' -T UnifiedGenotyper'
                               ' -R '+thisDir+'ref/ucsc.hg19.fasta'
                               ' -I '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.bam'
                               ' -o '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.unifiedGenotyper.vcf'
                               ' -dfrac 1'
                               ' -glm BOTH'
                               ' -minIndelFrac 0.01'
                               ' -mbq 10'
                               ' -L chr13:32889617-32973809'
                               ' -L chr17:41196312-41277500 -nt '+threads,
                               shell=True,stderr=sp.STDOUT)
    # Annotate variants
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal.unifiedGenotyper.ann.vcf') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal.unifiedGenotyper.ann.vcf.gz'):
        output=sp.check_output("java -jar "+configs[4]+"snpEff.jar hg19"
                               ' '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.unifiedGenotyper.vcf'
                               ' > '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.unifiedGenotyper.ann.vcf',
                               shell=True,stderr=sp.STDOUT)
    # Convert VCF to ANNOVAR input file
    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal.unifiedGenotyper.ann.avinput') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
                          '.sorted.read_groups.realigned.recal.unifiedGenotyper.ann.avinput.gz'):
        output=sp.check_output(configs[5]+"convert2annovar.pl"
                               ' -format vcf4 -includeinfo'
                               ' '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.unifiedGenotyper.ann.vcf'
                               ' >'+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.unifiedGenotyper.ann.avinput',
                               shell=True,stderr=sp.STDOUT)
##    # Recalculate allele distribution
##    if not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
##                          '.sorted.read_groups.realigned.recal.unifiedGenotyper.ann.ad_recal.avinput') and not os.path.exists(outDir+'patient_'+patNum+'/patient_'+patNum+''
##                          '.sorted.read_groups.realigned.recal.unifiedGenotyper.ann.ad_recal.avinput.gz'):
##        output=sp.check_output('python3 '+thisDir+'recalculateAlleleDistribution.py'
##                               ' -in '+outDir+'patient_'+patNum+'/patient_'+patNum+'.sorted.read_groups.realigned.recal.unifiedGenotyper.ann.avinput',
##                               shell=True,stderr=sp.STDOUT).decode('utf-8')
##        if 'ERROR' in output:
##            print(output)
##            exit(0)
    q1.put(1)
    showPercWork(q1.qsize(),allWork)

# Section of input arguments
par=argparse.ArgumentParser(description='This script do all processing of BRCA sequencing data')
par.add_argument('--readsFiles_r1','-r1',dest='readsFiles1',type=str,help='regular expression for choosing files with R1 reads',required=True)
par.add_argument('--readsFiles_r2','-r2',dest='readsFiles2',type=str,help='regular expression for choosing files with R2 reads (optional)',required=False)
par.add_argument('--readsFiles_N','-rN',dest='readsFilesN',type=str,help='regular expression for choosing files with native R1 reads (including Undetermined) to evaluate effectivity of trimming reads  (alternative). If you do not have trimmed reads, do not use this parameter',required=False)
par.add_argument('--patientsTable','-pat',dest='patientsTable',type=str,help='table with information about each patient: ngs_num patient_id barcode1 barcode2',required=True)
par.add_argument('--primersCoords','-primer',dest='primersCoords',type=str,help='table with information about amplicon coordinates without column headers: amplicon_number | chromosome | start | end. (Is not required)',required=False)
par.add_argument('--fixMisencodedQuals','-fix',dest='fixMisEncoded',action='store_true',help='this parameter is needed if GATK shows error of quality encoding')
par.add_argument('--outDir','-out',dest='outDir',type=str,help='directory for output',required=True)
par.add_argument('--threads','-th',dest='threads',type=int,help='number of threads',default=2)
par.add_argument('--tool-threads','-tt',dest='toolThreads',type=int,help='number of threads for each tool. Number of --threads multiplied by the number of --tool-threads must not exceed number of CPU cores',default=1)
par.add_argument('--run-name','-run',dest='runName',type=str,help='Name of run. This name will be added into the name of the final output file. Default: BRCA',required=False,default='BRCA')
par.add_argument('--without-joinment','-notjoin',dest='notToJoin',action='store_true',help='use this parameter if you only want to process patients reads separately without joining them (useful if you have an opportunity to separate processing onto several machines)')
par.add_argument('--only-join','-onlyjoin',dest='onlyJoin',action='store_true',help='use this parameter if you only want to join already processed patients reads')
par.add_argument('--language','-lang',dest='lang',type=str,help='Language of report and text on figures (russian or english). Default: english',default='english')
args=par.parse_args()

langs=['russian','english']
if args.lang not in langs:
    print('#'*10,'\nWARNING! Chosen language is not accepted. Use default english...')
if args.notToJoin and args.onlyJoin:
    print("ERROR: only one of parameters can be used, -notjoin or -onlyjoin. But you have tried to use both")
    exit(0)
readsFiles1=sorted(glob.glob(args.readsFiles1))
if len(readsFiles1)==0:
    print('ERROR: no files for reads R1 were selected!'); exit(0)
if args.readsFiles2:
    readsFiles2=sorted(glob.glob(args.readsFiles2))
    if len(readsFiles2)==0:
        print('ERROR: no files for reads R2 were selected!'); exit(0)
else:
    readsFiles2=[False]*len(readsFiles1)
if args.readsFilesN:
    readsFilesN=sorted(glob.glob(args.readsFilesN))
    if len(readsFilesN)==0:
        print('ERROR: no files for native reads were selected!'); exit(0)
##    if not (0<=len(readsFilesN)-len(readsFiles1)<=1):
##        print('ERROR: number of files with native R1 reads should be equal to or be one more than the number of files with trimmed R1 reads '); exit(0)
patientsTable=args.patientsTable
if not os.path.exists(patientsTable):
    print('ERROR: patients table file does not exist!'); exit(0)
if args.primersCoords and not os.path.exists(args.primersCoords):
    print('ERROR: primers coords table file does not exist!'); exit(0)
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
            exit(0)
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
            exit(0)

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
                    exit(0)
            patNums.append(m[0])
        # We have index numbers in our reads files names
        elif len(m)>0 and len(m2)>0:
            patNums.append(m2[0])
        else:
            print('INTERNAL ERROR: There is no numbers in the names of reads files!')
            exit(0)
else:
    patNums=['1']
allWork=len(patNums)
outDirs=[outDir]*allWork
if not args.onlyJoin:
    # So on this step we start from mapping reads
    q1=Queue()
    p=Pool(t)
    # Determine if there has been an error of quality encoding
    if args.fixMisEncoded:
        fixMis=[True]*len(readsFiles1)
    else:
        fixMis=[False]*len(readsFiles1)
    print('Process reads...')
    showPercWork(0,allWork)
    p.map(processPatient,zip(readsFiles1,readsFiles2,outDirs,patNums,fixMis))
    print('\nDone')
if args.notToJoin:
    exit(0)
elif args.onlyJoin:
    ds=glob.glob(outDir+'patient_*/*.unifiedGenotyper.ann.avinput')
    if len(ds)<len(patNums):
        print('ERROR: Number of files with completely processed reads is less than number of patients!')
        print('Number of patients:',len(patNums))
        print('Numner of files with completely processed reads (*.unifiedGenotyper.ann.avinput):',len(ds))
        print('Try to start analysis again without -onlyjoin')
        exit(0)
# Join mutations of all patients to one file
print('Joining mutations of all patients...')
output=sp.check_output("python3 "+thisDir+"joinMutations.py -v '"+outDir+"patient_*/*.unifiedGenotyper.ann.avinput' -o "+outDir+'allPatients.unifiedGenotyper.ann.avinput',shell=True,stderr=sp.STDOUT)
if 'ERROR' in str(output):
    print(output.decode('utf-8'))
    exit(0)
print('Done')
# Annotate them with ANNOVAR
print('Annotating mutations of all patients with ANNOVAR...')
output=sp.check_output(configs[5]+"table_annovar.pl "+outDir+'allPatients.unifiedGenotyper.ann.avinput '+configs[5]+'humandb/ -buildver hg19 -out '+outDir+'allPatients.unifiedGenotyper.ann.avinput.ann -remove -protocol refGene,cosmic70,esp6500siv2_all,exac03,kaviar_20150923,1000g2015aug_all,avsnp147,clinvar_20170905,ljb26_all -operation g,f,f,f,f,f,f,f,f -nastring . -otherinfo',shell=True,stderr=sp.STDOUT)
if 'ERROR' in str(output):
    print(output.decode('utf-8'))
    exit(0)
print('Done')
# Add positions in old references
print('Adding positions in old references...')
output=sp.check_output("python3 "+thisDir+"addPosOldRef.py -i "+outDir+'allPatients.unifiedGenotyper.ann.avinput.ann.hg19_multianno.txt -p '+patientsTable,shell=True,stderr=sp.STDOUT)
if 'ERROR' in str(output):
    print(output.decode('utf-8'))
    exit(0)
print('Done')
# Annotate by Clinical information
print('Annotating with BIC...')
output=sp.check_output("python3 "+thisDir+"annotateClinVars.py -i "+outDir+'allPatients.unifiedGenotyper.ann.avinput.ann.hg19_multianno.withOurCoordinates.xls',shell=True,stderr=sp.STDOUT)
if 'ERROR' in str(output):
    print(output.decode('utf-8'))
    exit(0)
print('Done')
# Convert to Excel-format
print('Converting to Excel-format...')
output=sp.check_output("python3 "+thisDir+"convertResultToExcel.py"
                       ' -i '+outDir+'allPatients.unifiedGenotyper.ann.avinput.ann.hg19_multianno.withOurCoordinates.clinical.xls'
                       ' -o '+outDir+'allPatients.'+args.runName+'.unifiedGenotyper.ann.avinput.ann.hg19_multianno.withOurCoordinates.clinical.excel.xls',
                       shell=True,stderr=sp.STDOUT)
if 'ERROR' in str(output):
    print(output.decode('utf-8'))
    exit(0)
print('Done')
if args.primersCoords:
    # Evaluating coverage
    print('Evaluating coverage...')
    process=sp.Popen(["python3",thisDir+"countCovForAmplic.py","min",
                      args.primersCoords,patientsTable,
                      ''+outDir+'patient_*/*.sorted.read_groups.realigned.recal.bam',
                      outDir+args.runName+'_amplicon_coverage_min.xls',str(t*int(args.toolThreads))],
                     stdout=sp.PIPE,shell=False,universal_newlines=True)
    for c in iter(process.stdout.readline,''):
        sys.stdout.write("\r"+str(c).replace('\n',''))
        sys.stdout.flush()
    print()
    process=sp.Popen(["python3",thisDir+"countCovForAmplic.py","mean",
                      args.primersCoords,patientsTable,
                      ''+outDir+'patient_*/*.sorted.read_groups.realigned.recal.bam',
                      outDir+args.runName+'_amplicon_coverage_mean.xls',str(t*int(args.toolThreads))],
                     stdout=sp.PIPE,shell=False,universal_newlines=True)
    for c in iter(process.stdout.readline,''):
        sys.stdout.write("\r"+str(c).replace('\n',''))
        sys.stdout.flush()
    print()
# Evaluate read statistics
print('Evaluating reads...')
if args.readsFilesN:
    output=sp.check_output('python3 '+thisDir+'getPercentOfProperlyTrimmedReads.py '
                           '-nat "'+args.readsFilesN+'" -trim "'+args.readsFiles1+'" -pat '+patientsTable+' -out '+outDir+args.runName+'.reads_statistics.xls',shell=True,stderr=sp.STDOUT)
else:
    output=sp.check_output('python3 '+thisDir+'getPercentOfProperlyTrimmedReads.py '
                           '-nat "'+args.readsFiles1+'" -pat '+patientsTable+' -out '+outDir+args.runName+'.reads_statistics.xls',shell=True,stderr=sp.STDOUT)
print('Done')
print('Evaluating uniformity of coverage...')
output=sp.check_output('python3 '+thisDir+'drawUniformityFigure.py '
                       '-cov '+outDir+args.runName+'_amplicon_coverage_mean.xls '
                       '-out '+outDir+args.runName+'_amplicon_coverage_mean.uniformity.tiff '
                       '-lang '+args.lang,shell=True,stderr=sp.STDOUT)
print('Done')
print('Creating report...')
output=sp.check_output('python3 '+thisDir+'makeReport.py '
                        '-read '+outDir+args.runName+'.reads_statistics.xls '
                       '-cov '+outDir+args.runName+'_amplicon_coverage_mean.xls '
                       '-res '+outDir+'allPatients.'+args.runName+'.unifiedGenotyper.ann.avinput.ann.hg19_multianno.withOurCoordinates.clinical.excel.xls '
                       '-fig '+outDir+args.runName+'_amplicon_coverage_mean.uniformity.tiff '
                       '-out '+outDir+args.runName+'.report.docx '
                       '-lang '+args.lang,shell=True,stderr=sp.STDOUT)
print('Done')

