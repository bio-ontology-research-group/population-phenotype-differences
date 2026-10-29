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
    plt.close(fig)
    paired=a.res/'paired_source_analysis/paired_summary.csv'
    if paired.exists():
        rows=list(csv.DictReader(paired.open()))
        data={(r['subset'],r['outcome']):r for r in rows if r['region']=='ALL'}
        fig,axes=plt.subplots(1,2,figsize=(7.1,2.25),gridspec_kw={'width_ratios':[1,1]})
        subsets=['all_publications','one_individual_each_output']
        for ax,outcome,scale,title,letter in zip(axes,['delta_n_terms','delta_test_share'],[1,100],
                ['Additional HPO terms','Change in test share (percentage points)'],['a','b']):
            for y,subset in enumerate(subsets[::-1]):
                r=data[subset,outcome];m,lo,hi=[float(r[k])*scale for k in ['mean','lo','hi']]
                ax.errorbar(m,y,xerr=[[m-lo],[hi-m]],fmt='o',ms=5,color='#0072B2',capsize=3,lw=1.2)
                ax.annotate('n = '+r['publications'],(m,y),xytext=(0,10),textcoords='offset points',ha='center',fontsize=7)
            ax.axvline(0,color='#888888',ls='--',lw=.7);ax.set_ylim(-.5,1.6)
            ax.set_xlabel(title);ax.text(-.04,1.02,letter,transform=ax.transAxes,fontweight='bold',fontsize=11)
        axes[0].set_yticks([0,1],['One individual\nin each output','All publications'])
        axes[1].set_yticks([0,1],[])
        fig.text(.54,.02,'Full text minus abstract; publication bootstrap 95% intervals. Clinical validation pending.',ha='center',fontsize=7)
        fig.subplots_adjust(left=.2,right=.98,bottom=.3,top=.86,wspace=.3)
        fig.savefig(out/'paired_source_effects.pdf',bbox_inches='tight');fig.savefig(out/'paired_source_effects.png',dpi=180,bbox_inches='tight')

if __name__=='__main__':main()
