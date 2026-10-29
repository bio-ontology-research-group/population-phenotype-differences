#!/usr/bin/env python3
"""Materialize explicitly reviewed binary tables; retain source-specific provenance.

Mappings below were checked against frozen HPO and the source text by the assistant.
This is not independent clinician validation. New tables require explicit curation.
"""
import argparse,csv,hashlib,json,xml.etree.ElementTree as ET
from pathlib import Path
from hpo_map2 import load_obo
from source_robustness import write

TABLES={'29482518':dict(gene='KMT2D',disease_id='OMIM:147920',table_id='Tab1',
    cohort_id='29482518_primary_two_boys',population_as_stated='Chinese',
    care_country='China',care_basis='our_institution: First Affiliated Hospital of Dalian Medical University',
    clinical_review='pending',source='https://doi.org/10.1186/s12881-018-0545-5',
    flags=['Case-report ascertainment; selected diagnostic features, two unrelated boys.',
           'Patient 2 age is five months in the table and six months in the text; retain both.',
           'Figure caption gives c.7948G>A, whereas abstract/table/body give c.7048G>A; retain variant discrepancy.',
           'Patient overlap with independent clinical evaluation requires checking.'],
    mappings={
        'Elongated palpebral fissures':('HP:0000637','equivalent'),
        'Eversion of the lateral third of the lower eyelid':('HP:0007655','equivalent'),
        'Laterally sparse eyebrows':('HP:0005338','equivalent'),
        'Depressed nasal tip':('HP:0000437','equivalent'),
        'Micrognathia':('HP:0000347','equivalent'),
        'Abnormal dentition':('HP:0000164','equivalent'),
        'Spinal column abnormalities':('HP:0000925','broad_source_term'),
        'Joint hypermobility/dislocation':('HP:0001382','Patient 2 text explicitly states joint hypermobility; negative table cell for Patient 1'),
        'Single palmar crease':('HP:0000954','equivalent'),
        'Clinodactyly of fifth digits':('HP:0030084','broader_mapping: table does not specify fingers versus toes'),
        'Postnatal growth deficiency':('HP:0008897','equivalent'),
        'Genitourinary anomalies':('HP:0000119','broad_source_term'),
        'Cerebellar vermis dysplasia':('HP:0002195','equivalent'),
        'Hypotonia':('HP:0001252','equivalent'),
        'Motor delay':('HP:0001270','equivalent'),
        'Feeding problem':('HP:0011968','equivalent'),
        'Recurrent infection':('HP:0002719','equivalent'),
    }),
    '25896430':dict(gene='KMT2D',disease_id='OMIM:147920',table_id='Tab1',
    cohort_id='25896430_author_causal_KMT2D_subset',population_as_stated='Han Chinese',
    care_country='China',care_basis='explicit: Beijing and Chongqing clinical records',
    clinical_review='pending',source='https://pubmed.ncbi.nlm.nih.gov/25896430/',
    source_patient_n=8,patient_columns=[1,3,4,6,8],
    flags=['Eight clinically diagnosed cases; only patients 1,3,4,6,8 assigned causal KMT2D variants by authors.',
           'Author missense classifications used computational evidence; retain distinction from modern clinical reclassification.',
           'Patients 4,6,8 carry reported truncating variants; a restricted sensitivity can use retained patient-level assessments.',
           'Phenotype-specific blank cells are unassessed; no whole-cohort denominator substitution.',
           'Patient overlap with independent clinical evaluation requires checking.'],
    mappings={
        'High/sparse eyebrows':('HP:0000534','broader_mapping: disjunctive source description'),
        'Long palpebral fissures':('HP:0000637','equivalent'),
        'Blue sclera':('HP:0000592','equivalent'),
        'Strabismus':('HP:0000486','equivalent'),
        'Ptosis':('HP:0000508','equivalent'),
        'Large ears':('HP:0000400','equivalent'),
        'Depressed nasal tip':('HP:0000437','equivalent'),
        'Micrognathia':('HP:0000347','equivalent'),
        'Abnormal dentition':('HP:0000164','equivalent'),
        'Clinodactyly of the 5th finger':('HP:0004209','equivalent'),
        'Hip dislocation':('HP:0002827','equivalent'),
        'Hyperlaxity':('HP:0001382','joint context in clinical description'),
        'Scoliosis/vertebral malformation':('HP:0000925','broader_mapping: scoliosis or vertebral malformation'),
        'Fingertip pads':('HP:0001212','equivalent'),
        'Mental retardation':('HP:0001249','historical source terminology'),
        'Hypotonia':('HP:0001252','equivalent'),
        'Seizures':('HP:0001250','equivalent'),
        'Microcephaly':('HP:0000252','equivalent'),
        'Postnatal growth retardation':('HP:0008897','equivalent'),
        'Feeding difficulties':('HP:0011968','equivalent'),
        'Frequent infections':('HP:0002719','equivalent'),
        'Cardiac anomalies':('HP:0001627','broad_source_term'),
    })}

def binary_counts(cells):
    known={'+' :1,'-':0,'−':0,'NA':None,'ND':None,'?':None,'':None}
    if any(c not in known for c in cells):raise ValueError('Unrecognized assessment symbol')
    vals=[known[c] for c in cells]
    return sum(v for v in vals if v is not None),sum(v is not None for v in vals),sum(v is None for v in vals)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--obo',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    T=load_obo(str(a.obo));rows=[];patient_rows=[];sources=[]
    for pmid,config in TABLES.items():
        source=a.raw/f'{pmid}.xml';root=ET.parse(source).getroot()
        table=next(t for t in root.findall('.//table-wrap') if t.get('id')==config['table_id'])
        source_n=config.get('source_patient_n',2)
        selected=config.get('patient_columns',list(range(1,source_n+1)))
        found=set()
        for tr in table.findall('.//tr'):
            cells=[' '.join(''.join(c.itertext()).split()) for c in tr if c.tag in ('th','td')]
            # Row-spanning category headings add one leading cell to some rows.
            labels=[i for i,c in enumerate(cells[:2]) if c in config['mappings']]
            if not labels:continue
            if len(labels)!=1:raise ValueError('Ambiguous phenotype row')
            offset=labels[0];label=cells[offset]
            vals=cells[offset+1:offset+1+source_n]
            if len(vals)!=source_n:raise ValueError('Unexpected number of patient columns')
            vals=[vals[i-1] for i in selected]
            found.add(label);hpo,note=config['mappings'][label]
            assert hpo in T
            k,n,u=binary_counts(vals)
            rows.append(dict(pmid=pmid,cohort_id=config['cohort_id'],gene=config['gene'],disease_id=config['disease_id'],
                hpo_id=hpo,hpo_label=T[hpo]['name'],affected=k,assessed=n,unassessed=u,cohort_n=len(selected),
                denominator_basis='explicit_binary_patient_table',source_location='Table 1: '+label,
                mapping_note=note,population_as_stated=config['population_as_stated'],care_country=config['care_country'],
                care_basis=config['care_basis'],curation_status='assistant_source_checked',clinical_review='pending',
                independent_evaluation_overlap_checked=False,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
            for i,status in zip(selected,vals):
                patient_rows.append(dict(pmid=pmid,cohort_id=config['cohort_id'],patient_label=f'Patient {i}',hpo_id=hpo,
                                         status='present' if status=='+' else 'absent' if status in ('-','−') else 'unassessed'))
        if found!=set(config['mappings']):raise ValueError('Missing curated table rows')
        sources.append(dict(pmid=pmid,**{k:v for k,v in config.items() if k!='mappings'}))
    write(a.out/'frequency_annotations.csv',rows);write(a.out/'patient_assessments.csv',patient_rows)
    (a.out/'source_review.json').write_text(json.dumps(sources,indent=2)+'\n')
    print(len(rows),'source-checked frequency rows; clinician review and overlap audit pending')

if __name__=='__main__':main()
