# This script creates figure with uniformity of coverage

import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import statistics as stat

enLang={'Coverage':'Coverage','MedAmplCov':'Median \nAmplicon Coverage',
        'AmplNum':'Amplicon Number','VarCoef':'Coefficient of Variation',
        'VarCoefCovAmpl':'Coefficient of Variation\nof Amplicon Coverage'}
ruLang={'Coverage':'Покрытие','MedAmplCov':'Медианное \nпокрытие ампликона',
        'AmplNum':'Номер ампликона','VarCoef':'Коэффициент вариации',
        'VarCoefCovAmpl':'Коэффициент вариации\nпокрытия ампликона'}

def set_style():
    plt.style.use(['seaborn-white', 'seaborn-paper'])
    matplotlib.rc("font", family="sans-serif",weight='bold')
    matplotlib.rc("text",color='black')

par=argparse.ArgumentParser(description='This script creates report about BRCA-analyzer results')
par.add_argument('--cov-stat-file','-cov',dest='covStatFile',type=str,help='file with statistics of coverage',required=True)
par.add_argument('--output-file','-out',dest='outFile',type=str,help='file for output',required=True)
par.add_argument('--language','-lang',dest='lang',type=str,help='language of report (russian or english). Default: english',default='english')
args=par.parse_args()

set_style()
langs=['russian','english']
if args.lang not in langs:
    print('#'*10,'\nWARNING! Chosen language is not accepted. Use default english...')
if args.lang=='russian': lang=ruLang
else: lang=enLang
file=open(args.covStatFile)
pats=[]
covList=[]
for string in file:
    if 'amplicon#' in string: continue
    cols=string.replace('\n','').split('\t')
    if 'DEL_' in cols[1] or 'empty' in cols[1] or cols[1]=='': continue
    pat=cols[0].replace('patient_','')
    pats.append(pat)
    covs0=list(map(float,cols[5:]))
    covs=[]
    for i,cov in enumerate(covs0):
        covs.append(cov)
    covList.append(covs)
data=np.array(covList)
maxi,maxj=data.shape
meanAmplCovs=[]
cvs=[]
for j in range(maxj):
    meanAmplCovs.append(round(float(stat.median(data[:,j])),1))
    cvs.append(round(float(stat.stdev(data[:,j])/stat.mean(data[:,j])),3))

x=list(range(1,len(meanAmplCovs)+1))
fig,ax1=plt.subplots(figsize=(15,4))
ax1.set_xlim(0,190)
pl1=ax1.bar(x,height=meanAmplCovs,label=lang['Coverage'])
ax1.set_ylim([0,None])
ax1.set_yticklabels(['{:3.0f}'.format(i) for i in ax1.get_yticks()],fontsize=16,fontweight='normal')
ax1.set_xticklabels(['{:3.0f}'.format(i) for i in ax1.get_xticks()],fontsize=16,fontweight='normal')
ax1.set_ylabel(lang['MedAmplCov'],fontsize=16,fontweight='bold')
ax1.set_xlabel(lang['AmplNum'],fontsize=16,fontweight='bold')
ax2=ax1.twinx()
pl2=ax2.plot(x,cvs,'r-',label=lang['VarCoef'])
y=[0,0.50,1.00,1.50,2.00,2.50]
ax2.set_yticklabels(['{:3.0f}%'.format(i*100) for i in y],fontsize=16,fontweight='normal')
ax2.set_xticklabels(['{:3.0f}'.format(i) for i in ax2.get_xticks()],fontsize=16,fontweight='normal')
ax2.set_ylim([0,2.6])
ax2.set_xlim(0,190)
ax2.set_ylabel(lang['VarCoefCovAmpl'],fontsize=16,fontweight='bold')
plt.title('Coverage Uniformity',fontsize=16,fontweight='bold')
ax1.legend(bbox_to_anchor=(0., 0.95, 1., .10),loc=3,fontsize=16)
ax2.legend(bbox_to_anchor=(0., 0.95, 1., .10),loc=4,fontsize=16)
plt.tight_layout()
plt.savefig(args.outFile)
plt.close()
