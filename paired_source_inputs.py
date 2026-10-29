#!/usr/bin/env python3
"""Prepare a deterministic paired-source sample and cache open-access full texts.

Eligibility: one Store patient and one abstract-derived patient for a publication-
gene unit in the existing validation data; one selected unit per publication.
Sampling precedes full-text access checks. Availability failures remain in the
manifest. This does not establish that the two patient mentions are the same person.
"""
import argparse,collections,concurrent.futures,csv,hashlib,json,random,re,time
import urllib.parse,urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

EPMC='https://www.ebi.ac.uk/europepmc/webservices/rest'


def get(url):
    for i in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'rare-disease-source-study/1.0'}),timeout=45) as r:
                return r.read()
        except Exception:
            if i==2: raise
            time.sleep(2*(i+1))


def text(el):
    return ' '.join(' '.join(el.itertext()).split()) if el is not None else ''


def body_text(body):
    parts=[]
    for el in body.iter():
        if el.tag in ('title','p'):
            parts.append(text(el))
        elif el.tag=='table-wrap':
            caption=text(el.find('caption'))
            rows=[' | '.join(text(c) for c in row if c.tag in ('td','th')) for row in el.iter('tr')]
            parts.append('TABLE: '+caption+'\n'+'\n'.join(rows))
    return '\n\n'.join(p for p in parts if p)


def fetch(row, raw):
    pmid=row['pmid'];out=dict(row)
    meta_path=raw/f'{pmid}.metadata.json';xml_path=raw/f'{pmid}.xml'
    try:
        if not meta_path.exists():
            query=urllib.parse.quote(f'EXT_ID:{pmid} AND SRC:MED')
            data=json.loads(get(f'{EPMC}/search?query={query}&format=json&resultType=core'))
            hits=data.get('resultList',{}).get('result',[])
            if not hits: return dict(out,status='metadata_not_found')
            meta_path.write_text(json.dumps(hits[0],indent=2)+'\n')
        meta=json.loads(meta_path.read_text());pmcid=meta.get('pmcid','')
        out.update(pmcid=pmcid,open_access=meta.get('isOpenAccess',''),title=meta.get('title',''))
        if not pmcid: return dict(out,status='no_pmc_fulltext')
        if not xml_path.exists(): xml_path.write_bytes(get(f'{EPMC}/{pmcid}/fullTextXML'))
        root=ET.fromstring(xml_path.read_bytes());body=root.find('.//body')
        if body is None: return dict(out,status='no_body')
        abstract='\n\n'.join(text(x) for x in root.findall('.//article-meta/abstract'))
        if not abstract: abstract=re.sub('<[^>]+>','',meta.get('abstractText',''))
        if not abstract: return dict(out,status='no_abstract')
        fulltext=abstract+'\n\n'+body_text(body)
        if len(fulltext)<len(abstract)+500: return dict(out,status='insufficient_body')
        # Keep full source intact; oversized text needs explicit chunked processing.
        for name,content in [('abstract',abstract),('fulltext',fulltext)]:
            (raw/f'{pmid}.{name}.txt').write_text(content+'\n')
            out[name+'_chars']=len(content)
            out[name+'_sha256']=hashlib.sha256(content.encode()).hexdigest()
        return dict(out,status='ready' if len(fulltext)<=90000 else 'requires_chunking')
    except Exception as e:
        return dict(out,status='fetch_failed',error=type(e).__name__+': '+str(e)[:180])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--per-region',type=int,default=30);ap.add_argument('--seed',type=int,default=20261029)
    a=ap.parse_args();a.raw.mkdir(parents=True,exist_ok=True);a.out.mkdir(parents=True,exist_ok=True)
    with (a.root/'results/v4/validation/abs_vs_store_units.csv').open() as f:
        units=[r for r in csv.DictReader(f) if r['level']=='patient' and r['region'] in ['EUR','MENA','SASIA','LATAM','AFR','EASIA']]
    # Input level spelling is intentionally validated rather than guessed silently.
    if not units: raise ValueError('No eligible patient-level units; inspect archived level values')
    rng=random.Random(a.seed);selected=[];seen=set();counts={}
    for pop in ['EUR','MENA','SASIA','LATAM','AFR','EASIA']:
        sub=sorted([r for r in units if r['region']==pop],key=lambda r:(r['pmid'],r['gene']))
        rng.shuffle(sub);unique=[]
        for r in sub:
            if r['pmid'] not in seen:
                seen.add(r['pmid']);unique.append(r)
        counts[pop]=len(unique)
        selected.extend(dict(pmid=r['pmid'],gene=r['gene'],region=pop) for r in unique[:a.per_region])
    (a.out/'sample.json').write_text(json.dumps({'seed':a.seed,'eligible_publications':counts,'selected':selected},indent=2)+'\n')
    rows=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(fetch,r,a.raw) for r in selected]
        for f in concurrent.futures.as_completed(futures):
            r=f.result();rows.append(r);print(r['pmid'],r['region'],r['status'],flush=True)
    rows.sort(key=lambda r:(r['region'],r['pmid']))
    fields=sorted({k for r in rows for k in r})
    with (a.out/'manifest.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    (a.out/'availability.json').write_text(json.dumps(dict(collections.Counter(r['status'] for r in rows)),indent=2)+'\n')

if __name__=='__main__':main()
