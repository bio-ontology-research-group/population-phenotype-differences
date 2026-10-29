#!/usr/bin/env python3
"""Create independent, blinded-to-model reference annotation packets.

Text source and geography may be recognizable. Model predictions, assigned region,
paired document identity and statistical results are withheld from reviewer packets.
Reviewer return files start blank; the program never manufactures annotations.
"""
import argparse,csv,hashlib,html,json,random,zipfile
from pathlib import Path

PROTOCOL='''# Independent source-text validation

Please review the supplied documents independently before discussing cases with the
other reviewer. Model predictions and assigned regional groups have been withheld.
The source type and geography may remain apparent from the clinical text; this is
not guaranteed source-type or geographical blinding.

This is phase 1: construct a reference directly from the source. Do not consult
model predictions. Each document names a target gene. Restrict annotation to primary
patients affected by that gene in this report. Exclude background/cited patients,
unaffected relatives, and laboratory models. Preserve patient labels and distinguish
individual patients from aggregates. If identity cannot be resolved, mark it uncertain.

Add one row per finding to findings.csv. State patient label, the source phrase,
present / explicitly_absent / uncertain, and the best HPO identifier and label.
An unmentioned finding is unknown and must not be entered as absent. Add the basis
of observation (history, examination, laboratory, imaging, physiological_test,
pathology, mixed or unclear), together with a short source locator/quote. Identify
findings at the specificity the source actually supports.

Complete one row per patient or aggregate in patients.csv. Record the diagnosis,
causal gene, variant, allelic state, ancestry and place of care only where supported.
Keep ancestry and care location separate. Use assessment_complete=yes only after the
whole supplied document has been reviewed; mark incomplete sources explicitly.

For any explicit frequency statement, add a row to frequencies.csv with the affected
count, number actually assessed for that finding, cohort definition, age range where
available, and supporting table/text location. Never assume every cohort member was
assessed because the paper gives a total cohort size. Flag historical cohorts and
possible patient overlap with other reports. Leave uncertain denominators blank.

After both independent reviews are returned, disagreements will be assembled in an
adjudication sheet. Phase 2 will compare the frozen reference with extracted findings
and HPO mappings, score errors and missed findings, and verify identity alignment of
abstract and full-text patients. Neither model agreement nor an empty return file
counts as completed clinician validation.

The source documents are supplied for internal research review. Do not redistribute
the packet as part of the public analysis repository. Return your completed CSV files
and document notes through the project team's usual channel.
'''


def csvfile(path,fields,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest',type=Path,required=True);ap.add_argument('--raw',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--public-summary',type=Path,required=True)
    ap.add_argument('--seed',type=int,default=20261031);a=ap.parse_args()
    if a.out.exists() and any(a.out.iterdir()):raise ValueError('Refusing to overwrite reviewer files')
    a.out.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader(a.manifest.open()));rng=random.Random(a.seed);sample=[]
    for pop in ['EUR','MENA','SASIA','LATAM','AFR','EASIA']:
        eligible=sorted([r for r in rows if r['status']=='ready' and r['region']==pop],key=lambda r:r['pmid'])
        rng.shuffle(eligible);sample.extend(eligible[:4])
    docs=[dict(r,source=s) for r in sample for s in ['abstract','fulltext']];rng.shuffle(docs)
    for i,d in enumerate(docs):d['document_id']=f'D{i+1:03d}'
    key=[{k:d[k] for k in ['document_id','pmid','gene','region','source']} for d in docs]
    csvfile(a.out/'PRIVATE_alignment_key.csv',list(key[0]),key)
    for reviewer in ['Paul','Marwa']:
        folder=a.out/reviewer;folder.mkdir();(folder/'documents').mkdir()
        order=docs.copy();rng.shuffle(order)
        index=['<!doctype html><meta charset="utf-8"><title>Clinical source review</title><h1>Clinical source review</h1><p>Read protocol.txt first. Predictions are withheld.</p><ol>']
        for d in order:
            source=a.raw/f"{d['pmid']}.{d['source']}.txt"
            (folder/'documents'/f"{d['document_id']}.txt").write_text('Document: '+d['document_id']+'\nTarget gene: '+d['gene']+'\n\n'+source.read_text())
            index.append(f'<li><a href="documents/{d["document_id"]}.txt">{d["document_id"]}</a> — target gene {html.escape(d["gene"])}</li>')
        (folder/'index.html').write_text('\n'.join(index+['</ol>']))
        (folder/'protocol.txt').write_text(PROTOCOL)
        csvfile(folder/'documents.csv',['document_id','target_gene','assessment_complete','source_incomplete','notes'],
                [dict(document_id=d['document_id'],target_gene=d['gene'],assessment_complete='',source_incomplete='',notes='') for d in order])
        csvfile(folder/'patients.csv',['document_id','patient_label','unit','number_of_patients','diagnosis','causal_gene','variant','allelic_state','ancestry','care_country','care_evidence','identity_evidence','notes'],[])
        csvfile(folder/'findings.csv',['document_id','patient_label','source_phrase','status','hpo_id','hpo_label','observation_basis','source_locator','notes'],[])
        csvfile(folder/'frequencies.csv',['document_id','phenotype','hpo_id','affected','assessed','cohort_definition','age_range','source_locator','overlap_notes','notes'],[])
        with zipfile.ZipFile(a.out/f'{reviewer}_source_validation.zip','w',zipfile.ZIP_DEFLATED) as z:
            for p in sorted(folder.rglob('*')):
                if p.is_file():z.write(p,p.relative_to(a.out))
    csvfile(a.out/'adjudication_template.csv',['document_id','patient_label','field','paul_assessment','marwa_assessment','adjudicated_value','rationale','adjudicator','date'],[])
    summary={'seed':a.seed,'reviewers':['Paul','Marwa'],'publications':len(sample),'documents_per_reviewer':len(docs),
             'publications_by_region':{p:sum(d['region']==p for d in sample) for p in sorted({d['region'] for d in sample})},
             'status':'prepared_not_reviewed','model_predictions_withheld':True,
             'phase':'independent reference annotation before prediction audit; pairing and statistical labels withheld',
             'source_manifest_sha256':hashlib.sha256(a.manifest.read_bytes()).hexdigest()}
    a.public_summary.write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
