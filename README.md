# BRCA-analyzer
**BRCA-analyzer** is an automatic **workflow** for an analysis of **BRCA1/2** genes **NGS** data. It has been developed and tested on the reads from MiSeq of more than 900 samples. All found pathogenic variations were confirmed by Sanger's sequencing.
## Docker use
BRCA-analyzer can be run as a Docker (www.docker.com) image. In this way you only need to install Docker (https://docs.docker.com/install/) and download image of BRCA-analyzer:
```
docker pull aakechin/brca-analyzer

```
After that, run downloaded image:
```
docker run -it --entrypoint 'bash' aakechin/brca-analyzer
```
You will be in the main image dorectory, where you can find folders with all necessary programs: brca_analyzer, annovar, snpEff etc. First of all, you need to unzip ANNOVAR databases and reference human genoome (they are very large, so it will take some time):
```
gunzip annovar/humandb/*.gz
gunzip brca_analyzer/ref/*.gz
```
Now, you can try to start example analysis or upload to the container (that you've obtained from the image) your FASTQ-files:
```
cd brca_analyzer/example 
python3 ../brca_analyzer.py -r1 'reads_trimmed/patient_*.r1.ad_trimmed.trimmed.qual_trimmed.fastq.gz' -r2 'reads_trimmed/patient_*.r2.ad_trimmed.trimmed.qual_trimmed.fastq.gz' -rN 'reads/*R1*' -out reads_trimmed_analysis/ -th 3 -tt 2 -run EXAMPLE -lang english
```
For uploading files to a container, you can use the following commands in other terminal:
```
docker ps -a
```
to know, what is your container's name. It will be written in the last column.
```
docker cp <file or directory that you want to copy to the container> <container name>/<folder where you want to put your files>
```
And then you can run analysis with command like above.

## Dependency
BRCA-analyzer needs the following external tools:
* **samtools** (version >= 1.2) (download it from https://github.com/samtools/samtools/releases);
* **bwa** (version >= 0.7.10) (download it from https://sourceforge.net/projects/bio-bwa/files/);
* **picard** (version >= 2.0.1) (download it from https://broadinstitute.github.io/picard/);
* **GenomeAnalysisToolKit** (version=3.6 or 3.7) (download it from https://software.broadinstitute.org/gatk/download/);
* **snpEff** (download it from http://snpeff.sourceforge.net/download.html);
* **Annovar** (download it from http://annovar.openbioinformatics.org/en/latest/user-guide/download/);

Also BRCA-analyzer needs **Python** (version >= 3) and several external modules “argparse” for reading input arguments, "numpy" for working with arrays, "python-docx" for creating report, "xlrd" for reading EXCEL-tables, "xlsxwriter" for creating EXCEL-tables, and "matplotlib" for drawing graphs. To install them, run the following command:
```
pip3 install argparse numpy python-docx xlrd biopython xlsxwriter matplotlib
```
 For using snpEff you will also need hg19 database. To download it use the following command:
```
java -jar <path-to-snpEff>/snpEff.jar download hg19
```
For Annovar you will need the following databases: refGene, cosmic70, esp6500siv2_all, exac03, kaviar_20150923, 1000g2015aug_all, avsnp150, clinvar_20180603, ljb26_all. To download them run for each database: 
```
<path-to-ANNOVAR>/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar <database name> humandb/
```
Also, you will need to upload and index hg19 reference genome in the directory ref/ of BRCA-analyzer. Download hg19 reference from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit. After that, go the directory with BRCA-analyzer and run tool twoBitToFa from UCSC (http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa). twoBitToFa is included into BRCA-analyzer, so it can be started by:
``` 
./twoBitToFa ref/hg19.2bit ref/ucsc.hg19.fasta 
```
Then, extract one more archive with human genome parts:
```
cd ref/ && tar -xf human_g1k_v37_chr13+17.fasta.tar.gz && cd ../
```
Index and create sequences dictionary for the reference genome:
``` 
bwa index ref/ucsc.hg19.fasta
java -jar <path-to-picard>/picard.jar CreateSequenceDictionary R=ref/ucsc.hg19.fasta O=ref/ucsc.hg19.dict
samtools faidx ref/ucsc.hg19.fasta
```
## Installation
To install BRCA-analyzer run `python install.py` with paths to external programs as input arguments (written paths are only examples, replace them with yours; 0 for bwa means that it was added to the PATH):
```
python install.py -bwa 0 -sam ~/tools/samtools/ -bcf ~/tools/bcftools/ -picard ~/Downloads/picard/dist/ -gatk ~/Downloads/GenomeAnalysis-3.7/ -snpeff ~/tools/snpEff/ -annovar ~/Downloads/ANNOVAR/
```
## Use
### Quick start guide
To run analysis, first, you can prepare patients table file (it is not necessary) with the following columns (you can use file from example directory as a template):
* Sample number – this number should correspond the number from file with reads. For example, if FASTQ-file has name “3_S3_L001_R1_001.fastq.gz”, the sample should has number 3.
* Sample ID – this ID will be added to the result file. It can be some IDs that your laboratory use for incoming samples.
* First index – number of first index used for barcoding.
* Second index – number of the second index used for barcoding.

After that you can run the analysis with the following commands. If you used amplicon-based library preparation, you need to cut primer sequences from reads with another our tool cutPrimers (https://github.com/aakechin/cutPrimers). In the example/ directory we've already added trimmed read sequences. And then run **BRCA-analyzer**
```
cd example 
python3 ../brca_analyzer.py -r1 'reads_trimmed/patient_*.r1.ad_trimmed.trimmed.qual_trimmed.fastq.gz' -r2 'reads_trimmed/patient_*.r2.ad_trimmed.trimmed.qual_trimmed.fastq.gz' -rN 'reads/*R1*' -out reads_trimmed_analysis/ -th 3 -tt 2 -run EXAMPLE -lang english
```
### Arguments
```
-h, --help            show this help message and exit
  --readsFiles_r1 READSFILES1, -r1 READSFILES1 regular expression for choosing files with R1 reads
  --readsFiles_r2 READSFILES2, -r2 READSFILES2 regular expression for choosing files with R2 reads (alternative)
  --readsFiles_N READSFILESN, -rN READSFILESN regular expression for choosing files with native R1 reads (including Undetermined) to evaluate effectivity of trimming reads (alternative). If you do not have trimmed reads, do not use this parameter
  --patientsTable PATIENTSTABLE, -pat PATIENTSTABLE table with information about each patient: ngs_num patient_id barcode1 barcode2
  --primersCoords PRIMERSCOORDS, -primer PRIMERSCOORDS table with information about amplicon coordinates without column headers: amplicon_number | chromosome | start | end. (Is not required)
  --fixMisencodedQuals, -fix this parameter is needed if GATK shows error of quality encoding
  --outDir OUTDIR, -out OUTDIR directory for output
  --threads THREADS, -th THREADS number of threads
  --tool-threads TOOLTHREADS, -tt TOOLTHREADS number of threads for each tool. Number of –threads multiplied by the number of --tool-threads must not exceed number of CPU cores
  --run-name RUNNAME, -run RUNNAME Name of run. This name will be added into the name of the final output file. Default: BRCA
  --without-joinment, -notjoin use this parameter if you only want to process patients reads separately without joining them (useful if you have an opportunity to separate processing onto several machines)
  --only-join, -onlyjoin use this parameter if you only want to join already processed patients reads
  --language LANG, -lang LANG Language of report and text on figures (russian or english). Default: english
```
### Opportunity of faster analysis
If you have several processors (e.g. on the server), each with several cores, you can run the analysis in two steps but several times faster. First, you need to separate all FASTQ-files onto several groups of 4 file pairs (R1 and R2 reads). It can be done with bash-script (if you don't know how, contact me, and I'll send you my script). After that for each group of files you can start the analysis on the different processors with parameter -notjoin that make BRCA-analyzer not to join result files of all samples into one. Second, you need to start analysis with the same parameters, but without -notjoin and with -onlyjoin (as -r1 and -r2 parameters you can choose any group of files). This process can be done on a single processor. But in this case we recommend to create output directory manually in advance.
## Citation
Kechin A, Khrapov E, Boyarskikh U, Kel A, Filipenko M. BRCA-analyzer: Automatic workflow for processing NGS reads of BRCA1 and BRCA2 genes. Comput Biol Chem. 2018 Dec;77:297-306. doi: 10.1016/j.compbiolchem.2018.10.012. PMID: 30408727
## Keywords
BRCA1, BRCA2, NGS data analysis, automatic workflow
