#!/usr/bin/env python3
"""Add source-checked, explicitly assessed cohorts to the binary-table curation.

Retains patient and family units. Ambiguous findings and source contradictions are
recorded rather than converted to absence. No model-only proposal is accepted.
"""
import argparse,csv,hashlib,json
from pathlib import Path
from hpo_map2 import load_obo
from source_robustness import write
from curate_frequency_tables import binary_counts

COHORTS={
 '27536553':dict(gene='MPV17',disease_id='OMIM:256810',cohort_id='27536553_four_children',population_as_stated='Korean',
  care_country='South Korea',care_basis='clinical cohort described by treating institution',families=['1','2','2','3'],
  variants=['c.197T>A;c.293C>T','c.451dupC','c.451dupC','c.280G>C;c.293C>T'],
  rows=[('Cholestasis','HP:0001396','++++','Table 1'),('Hepatic failure','HP:0001399','++++','Table 1 and outcomes'),
        ('Developmental delay','HP:0001263','++++','Table 1 and Case descriptions'),('Hypotonia','HP:0001252','++++','Table 1 and Case descriptions'),
        ('Recurrent vomiting','HP:0002013','++++','Table 1'),('Lactic acidosis','HP:0003128','++++','Table 1 and Case descriptions'),
        ('Elevated AST','HP:0031956','++++','Table 1: explicit normal range and four values'),
        ('Elevated ALT','HP:0031964','++++','Table 1: explicit normal range and four values'),
        ('Elevated GGT','HP:0030948','++++','Case descriptions: elevated in all patients'),
        ('Elevated total bilirubin','HP:0003573','++++','Case descriptions: elevated in all patients'),
        ('Elevated direct bilirubin','HP:0002908','++++','Case descriptions: elevated in all patients'),
        ('Short stature','HP:0004322','---+','Table 1 height SDs; text defines below -2 SD'),
        ('Underweight','HP:0004325','-+-+','Table 1 weight SDs; text defines below -2 SD'),
        ('Seizure','HP:0001250','----','Case descriptions: no patient developed seizures')],
  flags=['Patients 2 and 3 are siblings; four patients represent three families.',
         'Severe hepatocerebral ascertainment; observed frequencies are conditional on this cohort.',
         'Brain MRI: table says patient 4 not done, whereas prose says all examined. Delayed myelination frequency deferred.',
         'Hypoglycemia reported in patients 1-3; assessment in patient 4 unspecified, so no 3/4 frequency assigned.',
         'Biopsy illustrations concern patient 1; no cohort-wide histology frequency inferred.']),
 '31666973':dict(gene='RPGRIP1',disease_id='OMIM:613826',cohort_id='31666973_four_Japanese_LCA_cases',population_as_stated='Japanese',
  care_country='Japan',care_basis='explicit: National Center for Child Health and Development',families=['1','1','1','2'],
  variants=['c.3565_3571delCGAAGGC','c.3565_3571delCGAAGGC','c.3565_3571delCGAAGGC','c.1467+1G>T;c.2710+374_2895+78del'],
  rows=[('Nyctalopia','HP:0000662','++++','Table 1'),('Nystagmus','HP:0000639','++++','Table 1'),
        ('Keratoconus','HP:0000563','----','Table 1'),('Cataract','HP:0000518','----','Table 1'),
        ('Macular atrophy','HP:0007401','----','Table 1: macular atrophy or staphyloma absent'),
        ('Chorioretinal atrophy','HP:0000533','++++','Table 1'),
        ('Undetectable light-adapted ERG','HP:0030465','++++','Table 1: cone responses not detectable'),
        ('Abnormal dark-adapted ERG','HP:0030469','++++','Table 1: rod responses reduced or subnormal'),
        ('Retinal thinning on OCT','HP:0030329','++++','Abstract and Results: ONL thinning in all four cases')],
  flags=['Cases 1-3 are siblings, including monozygotic twins; four patients represent two families.',
         'These patients were previously reported; count this follow-up cohort only once.',
         'Assessment ages differ by test, as specified in Table 1; frequencies summarize the reported assessment, not one common age.',
         'Attenuated vessels and peripheral pigmentation have undefined plus/minus symbols; defer those frequencies.',
         'Oculo-digital sign explicitly absent in all four; no suitable exact term found in frozen HPO, retained as unmapped source evidence.'])}

# Source-reported causal variants used only to audit overlap, not to reclassify them.
TABLE_VARIANTS={
 '29482518':{'Patient 1':'c.5235delA','Patient 2':'c.7048G>A'},
 '25896430':{'Patient 1':'c.16295G>A','Patient 3':'c.8639T>C','Patient 4':'c.3095delT','Patient 6':'c.4395dupC','Patient 8':'c.4140T>A'}}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tables',type=Path,required=True);ap.add_argument('--raw',type=Path,required=True)
    ap.add_argument('--obo',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True);T=load_obo(str(a.obo))
    rows=list(csv.DictReader((a.tables/'frequency_annotations.csv').open()))
    patients=list(csv.DictReader((a.tables/'patient_assessments.csv').open()))
    sources=json.loads((a.tables/'source_review.json').read_text());metadata=[]
    for source in sources:
        for label,variants in TABLE_VARIANTS[source['pmid']].items():
            metadata.append(dict(pmid=source['pmid'],cohort_id=source['cohort_id'],gene=source['gene'],patient_label=label,
                                 family_id=source['pmid']+':'+label,variants_as_reported=variants))
    for pmid,c in COHORTS.items():
        source=a.raw/f'{pmid}.xml';sha=hashlib.sha256(source.read_bytes()).hexdigest()
        for label,h,pattern,location in c['rows']:
            assert h in T and len(pattern)==4
            k,n,u=binary_counts(list(pattern))
            rows.append(dict(pmid=pmid,cohort_id=c['cohort_id'],gene=c['gene'],disease_id=c['disease_id'],
                hpo_id=h,hpo_label=T[h]['name'],affected=k,assessed=n,unassessed=u,cohort_n=4,
                denominator_basis='explicit_assessment_table_or_cohort_statement',source_location=location+': '+label,
                mapping_note='source_checked; exact or conservative broader HPO mapping',population_as_stated=c['population_as_stated'],
                care_country=c['care_country'],care_basis=c['care_basis'],curation_status='assistant_source_checked',
                clinical_review='pending',independent_evaluation_overlap_checked=False,source_sha256=sha))
            for i,s in enumerate(pattern,1):
                patients.append(dict(pmid=pmid,cohort_id=c['cohort_id'],patient_label=f'Patient {i}',hpo_id=h,
                                     status='present' if s=='+' else 'absent'))
        for i,(family,variant) in enumerate(zip(c['families'],c['variants']),1):
            metadata.append(dict(pmid=pmid,cohort_id=c['cohort_id'],gene=c['gene'],patient_label=f'Patient {i}',
                                 family_id=pmid+':'+family,variants_as_reported=variant))
        sources.append(dict(pmid=pmid,**{k:v for k,v in c.items() if k not in ['rows','variants','families']}))
    write(a.out/'frequency_annotations.csv',rows);write(a.out/'patient_assessments.csv',patients)
    write(a.out/'patient_metadata.csv',metadata)
    (a.out/'source_review.json').write_text(json.dumps(sources,indent=2)+'\n')
    print(len(rows),'frequency rows;',len(patients),'patient assessments;',len(metadata),'patients')

if __name__=='__main__':main()
