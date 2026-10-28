#!/usr/bin/env python3
"""Fetch PubMed case-report abstracts with affiliations (expansion corpus beyond PMC open access).

  python3 pubmed_abstracts.py --query Q --years 1980-2026 --out abstracts.jsonl [--exclude pmids.txt]

Splits the query by publication year (and by month when a year exceeds 9,999 records) to stay
under the E-utilities esearch limit; efetch XML in batches of 200; resumable.
Each line: {pmid, year, journal, title, abstract, first_country, countries, mesh, pubtypes}.
"""
import argparse, json, os, re, sys, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET

EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def http(url, data=None):
    for i in range(5):
        try:
            return urllib.request.urlopen(url, data=data, timeout=120).read()
        except Exception:  # noqa: BLE001
            time.sleep(3 * (i + 1))
    raise RuntimeError("eutils failed")


def search(term):
    d = json.loads(http(f"{EU}/esearch.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "term": term, "retmax": 9999, "retmode": "json"})))
    time.sleep(0.35)
    return int(d["esearchresult"]["count"]), d["esearchresult"]["idlist"]


def all_pmids(q, y0, y1):
    out = set()
    for y in range(y0, y1 + 1):
        n, ids = search(f"({q}) AND {y}[dp]")
        if n > 9999:
            for m in range(1, 13):
                n2, ids2 = search(f"({q}) AND {y}/{m:02d}[dp]"); out |= set(ids2)
        else:
            out |= set(ids)
        print(f"  {y}: {n} (total {len(out)})", flush=True)
    return sorted(out)



CTRY = ['Saudi Arabia', 'United Arab Emirates', 'Qatar', 'Kuwait', 'Oman', 'Bahrain', 'Yemen', 'Jordan', 'Lebanon', 'Syria', 'Iraq', 'Iran', 'Turkey', 'Türkiye', 'Egypt', 'Morocco', 'Tunisia', 'Algeria', 'Libya', 'Sudan', 'Pakistan', 'Afghanistan', 'Palestine', 'India', 'Bangladesh', 'Sri Lanka', 'Nepal', 'China', 'Japan', 'Korea', 'Taiwan', 'Hong Kong', 'Singapore', 'Malaysia', 'Indonesia', 'Thailand', 'Vietnam', 'Viet Nam', 'Philippines', 'Australia', 'New Zealand', 'USA', 'United States', 'Canada', 'Mexico', 'Brazil', 'Argentina', 'Chile', 'Colombia', 'Peru', 'United Kingdom', 'UK', 'England', 'Scotland', 'Wales', 'Ireland', 'France', 'Germany', 'Netherlands', 'Belgium', 'Switzerland', 'Austria', 'Italy', 'Spain', 'Portugal', 'Sweden', 'Norway', 'Denmark', 'Finland', 'Iceland', 'Poland', 'Czech', 'Hungary', 'Greece', 'Israel', 'Russia', 'Ukraine', 'Romania', 'Bulgaria', 'Serbia', 'Croatia', 'Slovenia', 'Slovakia', 'Estonia', 'Latvia', 'Lithuania', 'South Africa', 'Nigeria', 'Kenya', 'Ethiopia', 'Ghana', 'Uganda', 'Tanzania', 'Cameroon', 'Cyprus', 'Luxembourg', 'Malta']
ALIAS = {'United States': 'USA', 'England': 'United Kingdom', 'Scotland': 'United Kingdom', 'Wales': 'United Kingdom', 'UK': 'United Kingdom',
         'Korea': 'South Korea', 'Hong Kong': 'China', 'Czech': 'Czechia', 'Türkiye': 'Turkey', 'Viet Nam': 'Vietnam'}
CRE = re.compile('|'.join(r'\b' + re.escape(c) + r'\b' for c in sorted(CTRY, key=len, reverse=True)))


def parse(xml):
    for art in ET.fromstring(xml).findall('.//PubmedArticle'):
        pmid = art.findtext('.//PMID')
        title = ''.join(art.find('.//ArticleTitle').itertext()) if art.find('.//ArticleTitle') is not None else ''
        ab = ' '.join(((t.get('Label') + ': ') if t.get('Label') else '') + ''.join(t.itertext()) for t in art.findall('.//Abstract/AbstractText'))
        year = art.findtext('.//PubDate/Year') or (art.findtext('.//PubDate/MedlineDate') or '')[:4]
        cs = []
        for au in art.findall('.//AuthorList/Author'):
            a = ' '.join((x.text or '') for x in au.findall('.//Affiliation'))
            if a:
                f = [ALIAS.get(x, x) for x in CRE.findall(a)]; cs.append(f[-1] if f else None)
        yield {"pmid": pmid, "year": year, "journal": art.findtext('.//Journal/ISOAbbreviation') or '', "title": title, "abstract": ab,
               "first_country": (cs[0] if cs else None), "countries": sorted({c for c in cs if c}),
               "mesh": [d.findtext('DescriptorName') for d in art.findall('.//MeshHeading')],
               "pubtypes": [p.text for p in art.findall('.//PublicationType')]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", required=True); ap.add_argument("--years", default="1980-2026")
    ap.add_argument("--out", required=True); ap.add_argument("--exclude")
    a = ap.parse_args()
    y0, y1 = map(int, a.years.split("-"))
    idfile = a.out + ".pmids"
    if os.path.exists(idfile): ids = open(idfile).read().split()
    else:
        ids = all_pmids(a.query, y0, y1); open(idfile, "w").write("\n".join(ids))
    excl = set(open(a.exclude).read().split()) if a.exclude else set()
    done = {json.loads(l)["pmid"] for l in open(a.out)} if os.path.exists(a.out) else set()
    todo = [p for p in ids if p not in excl and p not in done]
    print(f"pmids {len(ids)}, excluded {len(set(ids) & excl)}, to fetch {len(todo)}", flush=True)
    with open(a.out, "a") as f:
        for i in range(0, len(todo), 200):
            xml = http(f"{EU}/efetch.fcgi", urllib.parse.urlencode({"db": "pubmed", "id": ",".join(todo[i:i + 200]), "retmode": "xml"}).encode())
            for r in parse(xml): f.write(json.dumps(r) + "\n")
            f.flush(); time.sleep(0.35)
            if (i // 200) % 25 == 0: print(f"  fetched {i + 200}/{len(todo)}", flush=True)
    print("done")


if __name__ == "__main__":
    main()
