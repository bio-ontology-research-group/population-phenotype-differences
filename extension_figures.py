#!/usr/bin/env python3
"""Figures for the v5 source/genetics/frequency extension, from saved results."""
import argparse,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--res',type=Path,required=True);a=ap.parse_args()
    out=a.res/'figures';out.mkdir(exist_ok=True)
    rows=list(csv.DictReader((a.res/'source_robustness/paired_ranking.csv').open()))
    data={(r['corpus'],r['population']):r for r in rows if r['resampling_unit']=='pmid'}
    pops=['EUR','MENA','SASIA','LATAM','AFR','EASIA']
    names=['EUR','MENA','South Asia','Latin America','Sub-Saharan Africa','East Asia']
    sources=[('store','Store','#0072B2','o'),('pmcp','PMC-Patients','#009E73','s'),('pubmed_abs','Abstracts','#D55E00','D')]
    plt.rcParams.update({'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,(a1,a2)=plt.subplots(1,2,figsize=(7.1,4.2),gridspec_kw={'width_ratios':[1.1,1]})
    for i,pop in enumerate(pops):
        y=5-i
        for offset,(corpus,label,color,marker) in zip([.22,0,-.22],sources):
            r=data[corpus,pop];x,lo,hi=[float(r[k]) for k in ['geometric_rank_ratio','ratio_lo','ratio_hi']]
            a1.errorbar(x,y+offset,xerr=[[x-lo],[hi-x]],fmt=marker,ms=4,color=color,capsize=2,lw=1)
            x,lo,hi=[float(r[k]) for k in ['top10_change_pp','top10_change_lo','top10_change_hi']]
            a2.errorbar(x,y+offset,xerr=[[x-lo],[hi-x]],fmt=marker,ms=4,color=color,capsize=2,lw=1)
    a1.set_xscale('log',base=2);a1.set_xlim(.5,16);a1.set_xticks([.5,1,2,4,8,16],['0.5','1','2','4','8','16'])
    a1.axvline(1,color='#999999',ls='--',lw=.7);a2.axvline(0,color='#999999',ls='--',lw=.7)
    a1.set_yticks(range(6),names[::-1]);a2.set_yticks(range(6),[])
    a1.set_xlabel('Geometric mean rank ratio\n(bedside only / full profile)')
    a2.set_xlabel('Change in top-10 recovery\n(percentage points)')
    for ax,label in [(a1,'a'),(a2,'b')]:
        ax.set_ylim(-.6,5.6);ax.text(-.04,1.02,label,transform=ax.transAxes,fontweight='bold',fontsize=11)
    fig.legend(handles=[Line2D([],[],marker=m,color=c,ls='',label=l) for _,l,c,m in sources],loc='lower center',ncol=3,frameon=False)
    fig.subplots_adjust(left=.22,right=.98,bottom=.24,top=.94,wspace=.12)
    fig.savefig(out/'source_specific_ranking.pdf',bbox_inches='tight');fig.savefig(out/'source_specific_ranking.png',dpi=180,bbox_inches='tight')

if __name__=='__main__':main()
