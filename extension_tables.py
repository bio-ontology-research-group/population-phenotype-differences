#!/usr/bin/env python3
"""Generate manuscript extension tables directly from archived analysis outputs."""
import argparse
import csv
from pathlib import Path


def read(p):
    return list(csv.DictReader(p.open()))


def table(caption,label,columns,headers,rows):
    return '\n'.join([r'\clearpage',r'\begin{table}[p]',r'\centering',
                      r'\caption{'+caption+'}',r'\label{'+label+'}',r'\small',
                      r'\begin{tabular}{'+columns+'}',r'\toprule',
                      ' & '.join(headers)+r' \\',r'\midrule',
                      *[' & '.join(map(str,row))+r' \\' for row in rows],
                      r'\bottomrule',r'\end{tabular}',r'\end{table}',''])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();v=a.root/'results/v5';parts=[]
    names={'store':'Store','pmcp':'PMC-Patients','pubmed_abs':'Abstracts'}
    rows=[]
    for r in read(v/'source_robustness/paired_ranking.csv'):
        if r['resampling_unit']!='pmid':continue
        interval=f"{float(r['geometric_rank_ratio']):.2f} [{float(r['ratio_lo']):.2f}, {float(r['ratio_hi']):.2f}]"
        rows.append([names[r['corpus']],r['population'],r['n_cases'],interval,
                     f"{float(r['top10_all']):.1f}",f"{float(r['top10_bed']):.1f}"])
    parts.append(table('Bedside restriction within each corpus and region. Rank ratios compare bedside-only with full profiles for the same cases; intervals resample publications (2,000 draws). Top-ten recovery is in percent. Reporting-group intervals are released separately.',
                       'tab:sourceranks','llrlrr',['Source','Region','$n$','Rank ratio [95\\% interval]','Full','Bedside'],rows))
    rows=[]
    labels={'delta_n_terms':'All terms','delta_n_test':'Test-based terms','delta_n_bed':'Bedside terms','delta_test_share':'Test share (points)'}
    rs=read(v/'paired_source_analysis/paired_summary.csv')
    # Output names are checked rather than silently omitting an intended outcome.
    share=next(k for k in {r['outcome'] for r in rs} if 'share' in k and k.startswith('delta'))
    labels[share]='Test share (points)'
    for subset in ['all_publications','one_individual_each_output']:
        selected=[r for r in rs if r['subset']==subset and r['region']=='ALL' and r['outcome'] in labels]
        if not selected:raise ValueError('Missing paired subset '+subset)
        for r in selected:
            factor=100 if 'share' in r['outcome'] else 1
            rows.append(['All' if subset=='all_publications' else 'One individual',labels[r['outcome']],r['publications'],
                         f"{factor*float(r['mean']):.2f}",f"[{factor*float(r['lo']):.2f}, {factor*float(r['hi']):.2f}]"])
    parts.append(table('Common-pipeline source effects: full-text minus abstract model outputs. Count outcomes retain empty outputs; test-share pairs require classified findings in both views. Intervals use 5,000 publication bootstrap draws. The one-individual subset is based on output cardinality; independent clinical identity alignment is pending.',
                       'tab:pairedsource','llrrl',['Subset','Outcome','Pairs','Difference','95\\% interval'],rows))
    parts.append(table('Supported exploratory gene--region comparisons after matching, with at least ten cases and three publications per group. None of the tested composition or term outcomes passed correction across its full analysis-version family. Zero indicates insufficient support rather than an estimated null effect.',
                       'tab:geneextension','lr',['Matching beyond gene','Eligible comparisons'],
                       [['Corpus',342],['Corpus and author region',3],['Corpus and care country, stated origin',0],
                        ['Corpus, disease and allelic state',6],['Corpus, disease and inheritance',7],
                        ['Corpus, protein consequence and allelic state',4],['Corpus, exact variant and allelic state',0]]))
    ann=read(v/'frequency_curation_combined/frequency_annotations.csv');meta=read(v/'frequency_curation_combined/patient_metadata.csv');rows=[]
    for p in ['29482518','25896430','27536553','31666973']:
        aa=[r for r in ann if r['pmid']==p];mm=[r for r in meta if r['pmid']==p]
        rows.append([r'\href{https://pubmed.ncbi.nlm.nih.gov/'+p+'/}{'+p+'}',aa[0]['gene'],aa[0]['care_country'],len(mm),len({r['family_id'] for r in mm}),len(aa)])
    parts.append(table('Source-checked frequency cohorts. Rows retain phenotype-specific assessed denominators; the patient count is not automatically the denominator. PMID links identify the source papers. Kabuki cohorts use the authors\' causal KMT2D subset; the truncating-only sensitivity excludes two computationally assigned missense cases in PMID 25896430. Cohort and assessment details, source discrepancies and patient-level assessments are released with the annotations. Independent clinician validation is pending.',
                       'tab:frequencycohorts','lllrrr',['PMID','Gene','Care country','Patients','Families','Rows'],rows))
    m=read(v/'frequency_benchmark/metrics.csv');rows=[]
    arms={'unweighted':'Unweighted','existing':'Existing frequencies','updated':'Updated frequencies',
          'association_only':'Added associations, unweighted','association_and_frequency':'Added associations and frequencies'}
    for arm,label in arms.items():
        r=next(r for r in m if r['subset']=='all' and r['arm']==arm)
        rows.append([label,r['n'],f"{float(r['mrr']):.6f}",f"{float(r['median_rank']):.1f}",f"{100*float(r['top10']):.1f}"])
    parts.append(table('Exploratory frequency ranking across 4,891 genes. All arms share query phenotypes, candidate genes and gene-derived information content. The first three arms share phenotype support; the final two separately add 12 positive disease--phenotype associations. MRR is mean reciprocal midrank; top ten is a percentage. Only three publications contribute queries. No sensitivity arm produced top-ten recovery.',
                       'tab:frequencyperformance','lrrrr',['Reference','Queries','MRR','Median rank','Top ten'],rows))
    rs=read(v/'frequency_benchmark/query_ranks.csv');lookup={(r['id'],r['arm']):r for r in rs};rows=[]
    for r in rs:
        if r['arm']!='existing':continue
        u=lookup[r['id'],'updated']
        rows.append([r['id'].replace('PAVS:',''),r['gene'],'Saudi' if r['corpus']=='pavs_clin' else 'US',
                     f"{float(r['rank']):g}",f"{float(u['rank']):g}",'Yes' if r['stringent_identity']=='True' else 'No'])
    parts.append(table('Every eligible clinical query in the fixed-support comparison. Lower ranks are better. The stringent identity sensitivity retains all seven Saudi queries and US query Q0000078; the other four US records have consistent KMT2D gene targets but discordant within-gene patient linkage across source tables. Original phenotypes are preserved. Tied-rank bounds and all sensitivity ranks are released in the machine-readable outputs.',
                       'tab:frequencyqueries','lllrrl',['Query','Gene','Cohort','Existing','Updated','Stringent'],rows))
    rows=[]
    for r in read(v/'frequency_source_contrast/within_disorder_comparison.csv'):
        rows.append([r['phenotype'].capitalize(),r['care_country'],r['patient_affected']+'/'+r['patient_assessed'],
                     r['concordant_family_positive']+'/'+r['families_assessed'],
                     f"[{100*float(r['family_ci_low']):.1f}, {100*float(r['family_ci_high']):.1f}]"])
    parts.append(table('Within-disorder source comparison for RPGRIP1-associated LCA6. German source: PMID 30576320, patients 14--17, four unrelated families; Japanese source: PMID 31666973, four patients in two families. Patient counts use explicit assessments. Family statuses are concordant for these traits; intervals are exact 95\\% binomial intervals (percent) across families. The German cataract-positive patient was aged 47; the Japanese table gives ages 1--9 with later tests at 7--15 years. Study, population and care country are confounded, precluding a population-effect estimate. This source comparison was kept separate from the frozen ranking annotation set.',
                       'tab:frequencycontrast','llrrl',['Finding','Care country','Patients','Families','Family interval'],rows))
    a.out.write_text('\n'.join(parts));print(a.out)


if __name__=='__main__':main()
