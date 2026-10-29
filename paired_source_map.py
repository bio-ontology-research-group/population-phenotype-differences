#!/usr/bin/env python3
"""Snapshot successful paired extractions into one shared phrase-mapping vocabulary."""
import argparse,collections,hashlib,json
from pathlib import Path


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--extract',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    data=a.extract.read_bytes();records={}
    # A producer may currently be appending its last line; use only complete lines.
    lines=data.splitlines(keepends=True)
    for l in lines:
        if not l.endswith(b'\n'):continue
        d=json.loads(l)
        if d.get('ok'):records[d['key']]=d
    phrases={};mentions=[];qc=[]
    for key,d in sorted(records.items()):
        for i,p in enumerate(d['result']['patients']):
            for status,field in [('present','observed_phenotypes'),('absent','excluded_phenotypes')]:
                for j,f in enumerate(p[field]):
                    if not isinstance(f,dict) or not isinstance(f.get('text'),str) or (f.get('hpo_label') is not None and not isinstance(f.get('hpo_label'),str)):
                        qc.append(dict(key=key,patient_index=i,finding_index=j,status=status,error='invalid_phrase_schema'))
                        continue
                    text=f['text'].strip();label=(f.get('hpo_label') or '').strip()
                    if not text:
                        qc.append(dict(key=key,patient_index=i,finding_index=j,status=status,error='empty_phrase'));continue
                    k=hashlib.sha256(json.dumps([text,label],ensure_ascii=False).encode()).hexdigest()
                    phrases[k]=dict(key=k,text=text,label=label)
                    mentions.append(dict(document_key=key,pmid=d['pmid'],gene=d['gene'],source=d['source'],
                                         patient_index=i,patient_label=p.get('local_id'),unit=p.get('unit'),
                                         number_of_patients=p.get('number_of_patients'),status=status,phrase_key=k,
                                         evidence=f.get('evidence')))
    for name,rows in [('phrases',phrases.values()),('mentions',mentions),('qc',qc)]:
        (a.out/f'{name}.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    (a.out/'snapshot.json').write_text(json.dumps({'successful_documents':len(records),'phrases':len(phrases),
        'mentions':len(mentions),'invalid_mentions':len(qc),'document_keys':sorted(records),
        'extraction_snapshot_sha256':hashlib.sha256(data).hexdigest()},indent=2)+'\n')
    print(len(records),'documents;',len(phrases),'unique phrases;',len(qc),'invalid mentions')

if __name__=='__main__':main()
