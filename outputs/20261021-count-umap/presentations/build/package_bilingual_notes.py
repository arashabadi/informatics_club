from pathlib import Path
import json,re,zipfile,xml.etree.ElementTree as E

ROOT=Path(__file__).resolve().parents[1]
ORIGINAL=Path('/Users/ashi/github/informatics_club/sctalk/08_from-count-to-umap/08_From-Count-Matrix-to-UMAP_2026-10-21.pptx')
EDITED=Path('/Users/ashi/Library/Application Support/Codex/artifact-sessions/exports/latest/d723f72c6c59a72f8821035b8b50c62ff65d471dfde9f7d8ec48ea9c1a9620df/presentation.pptx')
DEST=ROOT/'build/candidate-bilingual-notes.pptx'
NS={'a':'http://schemas.openxmlformats.org/drawingml/2006/main','p':'http://schemas.openxmlformats.org/presentationml/2006/main'}
for prefix,uri in NS.items():E.register_namespace(prefix,uri)
E.register_namespace('r','http://schemas.openxmlformats.org/officeDocument/2006/relationships')
def tag(n):return '{'+NS['a']+'}'+n
persian=re.compile('[\u0600-\u06ff]')

def format_notes(data):
    root=E.fromstring(data)
    for para in root.findall('.//a:p',NS):
        text=''.join(t.text or '' for t in para.findall('.//a:t',NS))
        fa=bool(persian.search(text));heading=text.strip().endswith(':') and len(text)<100
        ppr=para.find('a:pPr',NS)
        if ppr is None:ppr=E.Element(tag('pPr'));para.insert(0,ppr)
        ppr.set('algn','r' if fa else 'l');ppr.set('rtl','1' if fa else '0')
        # Explicit formatting makes the bilingual sections readable in Notes.
        props=[]
        default=ppr.find('a:defRPr',NS)
        if default is None:default=E.SubElement(ppr,tag('defRPr'))
        props.append(default)
        for run in para.findall('a:r',NS):
            rp=run.find('a:rPr',NS)
            if rp is None:rp=E.Element(tag('rPr'));run.insert(0,rp)
            props.append(rp)
        for rp in props:
            rp.set('sz','1200');rp.set('b','1' if heading else '0');rp.set('lang','fa-IR' if fa else 'en-US')
            if 'rtl' in rp.attrib:del rp.attrib['rtl']
            for key,font in [('latin','Arial'),('cs','Arial Unicode MS')]:
                node=rp.find('a:'+key,NS)
                if node is None:node=E.SubElement(rp,tag(key))
                node.set('typeface',font)
    return E.tostring(root,encoding='utf-8',xml_declaration=True)

expected=json.loads((ROOT/'build/bilingual-notes.json').read_text())
with zipfile.ZipFile(ORIGINAL) as src,zipfile.ZipFile(EDITED) as edit,zipfile.ZipFile(DEST,'w') as out:
    for i in range(1,35):
        rel=f'ppt/slides/_rels/slide{i}.xml.rels'
        targets=lambda raw:sorted((x.get('Type'),x.get('Target')) for x in E.fromstring(raw))
        assert targets(src.read(rel))==targets(edit.read(rel)),f'Unexpected slide relationship target change: {i}'
        nr=f'ppt/notesSlides/_rels/notesSlide{i}.xml.rels'
        assert targets(src.read(nr))==targets(edit.read(nr)),f'Unexpected notes target change: {i}'
        note_root=E.fromstring(edit.read(f'ppt/notesSlides/notesSlide{i}.xml'))
        assert not any(any(k.startswith('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}') for k in x.attrib) for x in note_root.iter()),f'Unexpected relationship-dependent content in note {i}'
    for info in src.infolist():
        data=src.read(info.filename)
        if re.fullmatch(r'ppt/notesSlides/notesSlide\d+.xml',info.filename):data=format_notes(edit.read(info.filename))
        out.writestr(info,data)

with zipfile.ZipFile(ORIGINAL) as src,zipfile.ZipFile(DEST) as out:
    changed=[n for n in src.namelist() if src.read(n)!=out.read(n)]
    assert len(changed)==34 and all(re.fullmatch(r'ppt/notesSlides/notesSlide\d+.xml',n) for n in changed)
    for i,x in enumerate(expected,1):
        root=E.fromstring(out.read(f'ppt/notesSlides/notesSlide{i}.xml'))
        lines=[''.join(t.text or '' for t in p.findall('.//a:t',NS)) for p in root.findall('.//a:p',NS)]
        note='\n'.join(lines).strip()
        assert note==x['notes'],f'Note content mismatch: slide {i}'
        for p in root.findall('.//a:p',NS):
            t=''.join(x.text or '' for x in p.findall('.//a:t',NS))
            if persian.search(t):
                pp=p.find('a:pPr',NS);assert pp.get('algn')=='r' and pp.get('rtl')=='1'
    print(json.dumps({'candidate':str(DEST),'notes_updated':34,'non_note_package_parts_changed':0,'bilingual_scripts_and_explanations_verified':34,'persian_RTL':'pass'},indent=2))
