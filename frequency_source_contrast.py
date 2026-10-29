#!/usr/bin/env python3
"""Source-verified within-disorder frequency comparison, separate from ranking.

German RPGRIP1 source reviewed after the first frozen ranking evaluation. It is
used for provenance comparison, not retrospectively added to the evaluated update.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
from scipy.stats import beta
from source_robustness import write


def interval(k,n):
    return (0. if k==0 else float(beta.ppf(.025,k,n-k+1)),
            1. if k==n else float(beta.ppf(.975,k+1,n-k)))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    source=a.root/'data/v5_private/frequency_sources/30576320.xml'
    tree=ET.parse(source).getroot()
    table=next(t for t in tree.findall('.//table-wrap') if t.get('id')=='pone.0205380.t003')
    # The four unrelated molecularly assigned RPGRIP1 patients are IDs 14--17.
    patients=[]
    for tr in table.findall('.//tbody/tr'):
        cells=[' '.join(' '.join(c.itertext()).split()) for c in tr]
        if cells and cells[0] in ['14','15','16','17']:
            assert 'RPGRIP1' in cells[2] and cells[3]=='LCA'
            assert all(cells[i] in ['yes','no','n.a.'] for i in [6,7,8])
            patients.append(dict(pmid='30576320',patient_label=cells[0],gene='RPGRIP1',
                                 disease_id='OMIM:613826',age_gender=cells[1],variant_as_reported=cells[2],
                                 family_id='30576320:'+cells[0],nystagmus=cells[6],cataract=cells[7],strabismus=cells[8],
                                 population_as_stated='German origin',care_country='Germany',
                                 care_setting='University Eye Hospital Tuebingen',source_location='Table 3; genotype Table 1',
                                 source_sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    assert len(patients)==4
    write(a.out/'german_patients.csv',patients)
    japanese=list(csv.DictReader((a.root/'results/v5/frequency_curation_combined/patient_assessments.csv').open()))
    meta={(r['pmid'],r['patient_label']):r for r in csv.DictReader((a.root/'results/v5/frequency_curation_combined/patient_metadata.csv').open())}
    rows=[]
    for trait,h in [('nystagmus','HP:0000639'),('cataract','HP:0000518')]:
        for source_id in ['30576320','31666973']:
            if source_id=='30576320':
                obs=[(r['family_id'],r[trait]=='yes') for r in patients if r[trait] in ('yes','no')]
            else:
                obs=[(meta[r['pmid'],r['patient_label']]['family_id'],r['status']=='present') for r in japanese
                     if r['pmid']==source_id and r['hpo_id']==h and r['status'] in ('present','absent')]
            families={}
            for f,s in obs:families.setdefault(f,set()).add(s)
            # No discordant families here; avoid implying any-member frequency
            # generally estimates the same quantity as patient prevalence.
            assert all(len(states)==1 for states in families.values())
            k=sum(s for f,s in obs);n=len(obs)
            kf=sum(next(iter(s)) for s in families.values());nf=len(families)
            lo,hi=interval(kf,nf)
            rows.append(dict(pmid=source_id,disease_id='OMIM:613826',hpo_id=h,phenotype=trait,
                             patient_affected=k,patient_assessed=n,concordant_family_positive=kf,
                             families_assessed=nf,family_fraction=kf/nf,family_ci_low=lo,family_ci_high=hi,
                             population='German origin' if source_id=='30576320' else 'Japanese',
                             care_country='Germany' if source_id=='30576320' else 'Japan',
                             age_information='current age 17,26,36,47 years' if source_id=='30576320' else
                             'table age 1,1,7,9 years; later test ages 7-15; clinical-trait assessment time not separately resolved',
                             clinical_review='pending'))
    write(a.out/'within_disorder_comparison.csv',rows)
    summary=dict(disorder='OMIM:613826',gene='RPGRIP1',comparable_phenotypes=2,
                 n_patients=8,n_independent_families=6,
                 frequency_estimand='Patient affected/assessed counts, plus exact 95% binomial intervals across families with concordant observed status.',
                 result='Nystagmus: 4/4 patients in each source; cataract: 1/4 German and 0/4 Japanese patients. The cataract-positive German patient is aged 47.',
                 identifiability='Study, reported population and care country are perfectly confounded; ages and variant distributions differ. No adjusted population coefficient is identifiable.',
                 scope='Descriptive provenance contrast; not a population prevalence comparison. Family intervals reflect the sampled families and do not remove ascertainment bias.',
                 genetic_assignment='German Table 1 reports likely-pathogenic biallelic variants for all four RPGRIP1 patients; no new variant classification.',
                 benchmark_separation='This source was reviewed after the frozen four-cohort ranking run and was not added to that benchmark.',
                 source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
