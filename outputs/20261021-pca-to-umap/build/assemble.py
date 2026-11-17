"""Retain source slides 1–3 exactly and all slide-4 content except its table.

Authored continuation uses the same source deck through Artifact Tool. Original
paper-figure slides are retained in an appendix. Relationship closures preserve
the source masters, logos, media, notes and editable objects.
"""
from pathlib import Path,PurePosixPath
import zipfile,xml.etree.ElementTree as E,posixpath,re,json,copy
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'build/source-preserved.pptx';AUTHORED=ROOT/'build/authored.pptx';OUT=ROOT/'build/candidate.pptx'
A='http://schemas.openxmlformats.org/drawingml/2006/main';P='http://schemas.openxmlformats.org/presentationml/2006/main';R='http://schemas.openxmlformats.org/officeDocument/2006/relationships';PKG='http://schemas.openxmlformats.org/package/2006/relationships';CT='http://schemas.openxmlformats.org/package/2006/content-types'
for pre,ns in [('a',A),('p',P),('r',R)]:E.register_namespace(pre,ns)
ns={'a':A,'p':P,'r':R}
def xml(node):return E.tostring(node,encoding='utf-8',xml_declaration=True)
def relpath(part):return str(PurePosixPath(part).parent/'_rels'/(PurePosixPath(part).name+'.rels'))
def resolve(part,target):return posixpath.normpath(target.lstrip('/') if target.startswith('/') else posixpath.join(posixpath.dirname(part),target))
src=zipfile.ZipFile(SOURCE);auth=zipfile.ZipFile(AUTHORED)
files={n:auth.read(n) for n in auth.namelist()}
types=E.fromstring(files['[Content_Types].xml']);st=E.fromstring(src.read('[Content_Types].xml'))
source_type={x.get('PartName').lstrip('/'):x.get('ContentType') for x in st if x.tag.endswith('Override')}
overrides={x.get('PartName'):x for x in types if x.tag.endswith('Override')}
defaults={x.get('Extension'):x.get('ContentType') for x in types if x.tag.endswith('Default')}
for t in st:
 if t.tag.endswith('Default') and t.get('Extension') not in defaults:
  types.append(copy.deepcopy(t));defaults[t.get('Extension')]=t.get('ContentType')
def settype(dest,original):
 if original not in source_type:return
 if '/'+dest in overrides:overrides['/'+dest].set('ContentType',source_type[original])
 else:
  el=E.SubElement(types,'{'+CT+'}Override',{'PartName':'/'+dest,'ContentType':source_type[original]});overrides['/'+dest]=el

appendix=[25,5,7,8,18,26,27,28,29,30,31,32,33,34,35]
source_slides={1:1,2:2,3:3,4:4,17:23,**{old:25+i for i,old in enumerate(appendix)}}
mapping={f'ppt/slides/slide{old}.xml':f'ppt/slides/slide{new}.xml' for old,new in source_slides.items()}
done=set()
def clone(part):
 dest=mapping.setdefault(part,'ppt/preserved/'+part.removeprefix('ppt/'))
 if part in done:return dest
 done.add(part);files[dest]=src.read(part);settype(dest,part)
 rp=relpath(part)
 if rp in src.namelist():
  rr=E.fromstring(src.read(rp))
  for rel in rr:
   if rel.get('TargetMode')=='External':continue
   target=resolve(part,rel.get('Target'));new=clone(target);rel.set('Target','/'+new)
  files[relpath(dest)]=xml(rr)
 return dest
for part in list(mapping):clone(part)

# Only the table XML on slide 4 is replaced. All other bytes remain identical.
original=src.read('ppt/slides/slide4.xml');table=E.fromstring(auth.read('ppt/slides/slide4.xml')).find('.//a:tbl',ns)
assert table is not None and len(table.findall('a:tr',ns))==7
# Keep the table within its original 1,152-pixel frame after adding a column.
grid=table.find('a:tblGrid',ns);widths=[.7,1,1,1,1,1.3];total_width=10972800
parts=[round(total_width*w/sum(widths)) for w in widths];parts[-1]+=total_width-sum(parts)
for col,width in zip(grid,parts):col.set('w',str(width))
tablebytes=E.tostring(table,encoding='utf-8')
match=re.search(rb'<a:tbl(?:\s[^>]*)?>.*?</a:tbl>',original,re.S);assert match
new4=original[:match.start()]+tablebytes+original[match.end():]
files['ppt/slides/slide4.xml']=new4
assert new4[:match.start()]==original[:match.start()]
assert new4[match.start()+len(tablebytes):]==original[match.end():]

# Renumber only the retained source slides beyond the protected opening.
for old,new in source_slides.items():
 if new<=4:continue
 part=f'ppt/slides/slide{new}.xml';rt=E.fromstring(files[part])
 for t in rt.findall('.//a:t',ns):
  if re.fullmatch(r'\d{2}(?: / Appendix)?',t.text or ''):t.text=f'{new:02}'+(' / Appendix' if new>=25 else '')
 files[part]=xml(rt)

pres=E.fromstring(files['ppt/presentation.xml']);rels=E.fromstring(files['ppt/_rels/presentation.xml.rels']);slist=pres.find('p:sldIdLst',ns)
maxid=max(int(x.get('id')) for x in slist);existingids={x.get('Id') for x in rels}
counter=10000
def rid():
 global counter
 while 'rId'+str(counter) in existingids:counter+=1
 s='rId'+str(counter);counter+=1;existingids.add(s);return s
for n in range(25,40):
 rr=rid();E.SubElement(rels,'{'+PKG+'}Relationship',{'Id':rr,'Type':R+'/slide','Target':f'/ppt/slides/slide{n}.xml'})
 maxid+=1;E.SubElement(slist,'{'+P+'}sldId',{'id':str(maxid),'{'+R+'}id':rr})
masters=pres.find('p:sldMasterIdLst',ns);maxmaster=max(int(x.get('id')) for x in masters)
for part in done:
 if source_type.get(part,'').endswith('slideMaster+xml'):
  rr=rid();E.SubElement(rels,'{'+PKG+'}Relationship',{'Id':rr,'Type':R+'/slideMaster','Target':'/'+mapping[part]})
  maxmaster+=1;E.SubElement(masters,'{'+P+'}sldMasterId',{'id':str(maxmaster),'{'+R+'}id':rr})
files['ppt/presentation.xml']=xml(pres);files['ppt/_rels/presentation.xml.rels']=xml(rels)
E.register_namespace('',CT)
files['[Content_Types].xml']=xml(types)

# Notes on protected slides remain exact. Direction is local in all new notes.
protected_notes=set()
for i in range(1,5):
 rr=E.fromstring(files[relpath(f'ppt/slides/slide{i}.xml')])
 for rel in rr:
  if rel.get('Type').endswith('/notesSlide'):protected_notes.add(rel.get('Target').lstrip('/'))
for part,data in list(files.items()):
 if 'notesSlide' not in part or not part.endswith('.xml') or '/_rels/' in part or part in protected_notes:continue
 rt=E.fromstring(data)
 for para in rt.findall('.//a:p',ns):
  text=''.join(t.text or '' for t in para.findall('.//a:t',ns));fa=bool(re.search('[\u0600-\u06ff]',text))
  pp=para.find('a:pPr',ns)
  if pp is None:pp=E.Element('{'+A+'}pPr');para.insert(0,pp)
  pp.set('algn','r' if fa else 'l');pp.set('rtl','1' if fa else '0')
  for rp in para.findall('a:r/a:rPr',ns):rp.set('lang','fa-IR' if fa else 'en-US')
 files[part]=xml(rt)

with zipfile.ZipFile(OUT,'w',compression=zipfile.ZIP_DEFLATED) as out:
 for n,data in files.items():out.writestr(n,data)
for i in [1,2,3]:assert files[f'ppt/slides/slide{i}.xml']==src.read(f'ppt/slides/slide{i}.xml')
report={'slides':39,'source_slides_1_3_byte_identical':True,'slide4_only_table_xml_changed':True,'core_route':'1–24','appendix':'25–39','cloned_source_parts':len(done),'source_slide_mapping':source_slides}
(ROOT/'build/preservation-report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
