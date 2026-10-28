"""Tag Phenopacket Store source PMIDs with first-author country from PubMed.
Usage: python3 -I tag_store_pmids_country.py store_pmids.txt out.tsv
Country detection is a regex over AffiliationInfo; see README.md for caveats."""
import sys,re,time,collections,urllib.request
import xml.etree.ElementTree as ET
pmids=sorted(set(open(sys.argv[1]).read().split()))
countries=['Saudi Arabia','United Arab Emirates','Qatar','Kuwait','Oman','Bahrain','Yemen','Jordan','Lebanon','Syria','Iraq','Iran','Turkey','Egypt','Morocco','Tunisia','Algeria','Libya','Sudan','Pakistan','India','Bangladesh','Sri Lanka','China','Japan','Korea','Taiwan','Hong Kong','Singapore','Malaysia','Indonesia','Thailand','Vietnam','Philippines','Australia','New Zealand','USA','United States','Canada','Mexico','Brazil','Argentina','Chile','Colombia','Peru','United Kingdom','UK','England','Scotland','Wales','Ireland','France','Germany','Netherlands','Belgium','Switzerland','Austria','Italy','Spain','Portugal','Sweden','Norway','Denmark','Finland','Iceland','Poland','Czech','Hungary','Greece','Israel','Russia','Ukraine','Romania','Bulgaria','Serbia','Croatia','Slovenia','Slovakia','Estonia','Latvia','Lithuania','South Africa','Nigeria','Kenya','Ethiopia','Ghana','Uganda','Tanzania','Cameroon']
alias={'United States':'USA','England':'United Kingdom','Scotland':'United Kingdom','Wales':'United Kingdom','UK':'United Kingdom','Korea':'South Korea','Hong Kong':'China','Czech':'Czechia'}
cre=re.compile('|'.join(r'\b'+re.escape(c)+r'\b' for c in sorted(countries,key=len,reverse=True)))
out=open(sys.argv[2],'w')
out.write('pmid\tyear\tjournal\tfirst_country\tall_countries\tn_authors_with_aff\tconsang_abstract\tmesh_geo\ttitle\n')
for i in range(0,len(pmids),150):
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id='+','.join(pmids[i:i+150])
    for attempt in range(4):
        try: xml=urllib.request.urlopen(url,timeout=120).read(); break
        except Exception: time.sleep(3*(attempt+1))
    for art in ET.fromstring(xml).findall('.//PubmedArticle'):
        pmid=art.findtext('.//PMID'); title=(art.findtext('.//ArticleTitle') or '').replace('\t',' ').replace('\n',' ')
        year=art.findtext('.//PubDate/Year') or art.findtext('.//PubDate/MedlineDate') or ''
        journal=art.findtext('.//Journal/ISOAbbreviation') or ''
        abstract=' '.join(t.text or '' for t in art.findall('.//AbstractText'))
        cons='yes' if re.search(r'consanguin',abstract,re.I) else 'no'
        cs=[]
        for au in art.findall('.//Author'):
            a=' '.join((x.text or '') for x in au.findall('.//Affiliation'))
            if a:
                f=[alias.get(x,x) for x in cre.findall(a)]; cs.append(f[-1] if f else None)
        allc=collections.Counter(c for c in cs if c)
        mesh=[d.findtext('DescriptorName') for d in art.findall('.//MeshHeading') if d.findtext('DescriptorName')]
        out.write('\t'.join([pmid,str(year)[:4],journal,(cs[0] if cs else None) or '',';'.join(f'{k}:{v}' for k,v in allc.most_common()),str(len(cs)),cons,';'.join(m for m in mesh if cre.search(m)),title])+'\n')
    time.sleep(0.4)
