#!/usr/bin/env python3
"""Audit frozen clinical queries against source diagnosis status before ranking.

Original patient phenotypes stay frozen. Gene-target identity is distinct from
unresolved within-gene patient/variant linkage in one published supplement.
"""
import argparse,collections,csv,hashlib,json,statistics
from pathlib import Path
from source_robustness import write

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pavs',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    normalized=a.pavs/'combined_normalized-gpt4oss.tsv'
    original={r['pavs_id']:r for r in csv.DictReader(normalized.open(),delimiter='\t')}
    index={r['id']:r for r in csv.DictReader((a.root/'results/phenopackets/index.tsv').open(),delimiter='\t')}
    frozen=list(csv.DictReader((a.root/'results/v5/frequency_sources/frozen_clinical_evaluation.csv').open()))
    rows=[]
    for r in frozen:
        s=original[r['id']];status=s['result_status'];eligible=status=='Solved'
        reason='reported_solved_source_record' if eligible else 'source_'+status.lower()
        subtype={'KMT2D':'OMIM:147920','RPGRIP1':'OMIM:613826','MPV17':'OMIM:256810'}[r['gene']]
        if r['id']=='PAVS:P0000527':
            eligible=False;reason='isolated_neuropathy_p.Arg41Gln_outside_hepatocerebral_update'
            subtype='OMIM:618400'
        linkage='source_record'
        if r['corpus']=='pavs_us':
            linkage='gene_target_consistent_across_supplement_tables; within_gene_patient_variant_linkage_uncertain'
            if r['id']=='PAVS:Q0000078':linkage='age_gene_and_id_concordant_between_source_tables'
        rows.append(dict(**r,source_id=s['source_id'],source_result_status=status,
                         source_diagnosis=s['disease_mondo_label'],source_variant_cdna=s['variant_cdna'],
                         source_variant_protein=s['variant_protein'],eligible_primary=eligible,eligibility_reason=reason,
                         benchmark_disorder=subtype,identity_check=linkage,
                         stringent_identity_sensitivity=eligible and (r['corpus']!='pavs_us' or r['id']=='PAVS:Q0000078'),
                         reference_overlap_status='pending'))
    write(a.out/'query_eligibility.csv',rows)
    counts=collections.Counter()
    for id,r in index.items():
        if r['corpus'] in ['pavs_clin','pavs_us']:
            counts[r['corpus'],original.get(id,{}).get('result_status','missing')]+=1
    write(a.out/'clinical_source_status.csv',[dict(corpus=c,result_status=s,n_cases=n) for (c,s),n in sorted(counts.items())])
    # The earlier nearest-case benchmark used clinical records with any gene label.
    # Its reference library excludes these routine-clinical corpora, so restricting
    # clinical query eligibility does not require changing the saved query ranks.
    rank_path=a.root/'results/v4/recovery/recovery_MENA_recovery_percase.csv'
    ranks=list(csv.DictReader(rank_path.open()));recovery=[]
    for corpus in ['pavs_clin','pavs_us']:
        cs=[r for r in ranks if r['corpus']==corpus]
        for subset in ['all_gene_labelled_queries','source_reported_solved']:
            sub=[r for r in cs if subset=='all_gene_labelled_queries' or original[r['id']]['result_status']=='Solved']
            vals=[float(r['rank']) for r in sub]
            recovery.append(dict(corpus=corpus,subset=subset,n_queries=len(vals),n_publications=len({r['pmid'] for r in sub}),
                                 top1=sum(x<=1 for x in vals)/len(vals),top10=sum(x<=10 for x in vals)/len(vals),
                                 mrr=sum(1/x for x in vals)/len(vals),median_rank=statistics.median(vals)))
    write(a.out/'existing_clinical_recovery.csv',recovery)
    summary=dict(frozen_candidates=len(rows),primary_eligible=sum(r['eligible_primary'] for r in rows),
                 stringent_identity_eligible=sum(r['stringent_identity_sensitivity'] for r in rows),
                 eligibility_reasons=dict(collections.Counter(r['eligibility_reason'] for r in rows)),
                 source_sha256=hashlib.sha256(normalized.read_bytes()).hexdigest(),
                 existing_recovery_sha256=hashlib.sha256(rank_path.read_bytes()).hexdigest(),
                 interpretation='Eligibility frozen before frequency ranking. Source-reported solved status is retained; this is not a new clinical variant classification. The US supplement links five phenotype records to five KMT2D diagnoses by ID; ages and the prose case illustrate within-gene identity discrepancies. Gene targets agree, and an identity-restricted sensitivity is specified.',
                 independence='Publication, reporting-group, variant and cohort-overlap checks required before ranking.')
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
