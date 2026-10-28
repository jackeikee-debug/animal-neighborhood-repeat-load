"""Offline evidence inventory checker; stdlib only, never executes analysis.

Default mode reads root manifests/directory metadata only. Content mode is
optional and explicit; this study runs it only on NEW tiny synthetic fixtures
and the NEW supplement package, not the already accepted scientific bundles.
Historical source_provenance and source-map paths are never followed.
"""
from pathlib import Path,PurePosixPath
from zipfile import ZipFile
import argparse,hashlib,json,re
CAP_MANIFEST=2097152;CAP_EXPANDED=268435456
def need(v,why):
 if not v:raise ValueError(why)
def safe_relative(value):
 need(isinstance(value,str)and value and '\\'not in value and '\x00'not in value,'invalid relativepath')
 need(not value.startswith('/')and not re.match(r'^[A-Za-z]:',value),'absolute/drivepath forbidden')
 need(all(x not in ['', '.', '..']for x in value.split('/')),'unsafe pathcomponent')
 need(not PurePosixPath(value).is_absolute(),'absolute pathforbidden');return value
def check_dependencies(archives):
 labels=[x['label']for x in archives];need(len(labels)==len(set(labels)),'duplicate archive label');known=set(labels);edges={x['label']:x.get('requires',[])for x in archives}
 need(all(r in known for v in edges.values()for r in v),'missing dependency')
 def visit(k,active,done):
  need(k not in active,'dependency cycle')
  if k in done:return
  active.add(k)
  for n in edges[k]:visit(n,active,done)
  active.remove(k);done.add(k)
 done=set()
 for k in labels:visit(k,set(),done)
def inspect_archive(path,manifest_sha256,mode='manifest',expected_manifest_bytes=None):
 need(mode in ['manifest','content'],'explicit mode')
 need(bool(re.fullmatch('[a-f0-9]{64}',manifest_sha256)),'invalid rootmanifestSHA')
 with ZipFile(path)as z:
  names=z.namelist();need(len(names)<=2000 and len(names)==len(set(names)),'duplicate/excessarchive members')
  candidates=[n for n in names if n.lower()=='manifest.json'];need(len(candidates)==1,'one rootmanifest required');root=candidates[0];info=z.getinfo(root)
  need(info.file_size<=CAP_MANIFEST,'manifest sizecap')
  if expected_manifest_bytes is not None:need(info.file_size==expected_manifest_bytes,'boundmanifestbytes')
  with z.open(root)as f:raw=f.read(CAP_MANIFEST+1)
  need(len(raw)<=CAP_MANIFEST and hashlib.sha256(raw).hexdigest()==manifest_sha256,'rootmanifestidentity')
  m=json.loads(raw);rows=m['files']if isinstance(m,dict)else m;need(isinstance(rows,list)and len(rows)<=1999,'manifest schema')
  paths=[];expanded=0
  for row in rows:
   p=safe_relative(row['path']);paths.append(p);n=row['bytes'];need(type(n)is int and n>=0,'nonnegative declaredbytes');expanded+=n;need(expanded<=CAP_EXPANDED,'expanded-bytecap')
   need(bool(re.fullmatch('[a-f0-9]{64}',row['sha256'])),'payloadSHA schema');need(p in names and z.getinfo(p).file_size==n,'missing member/bytesmismatch')
  sidecar_verified=False;allowed_extra=set()
  if 'MANIFEST.sha256'in names and 'MANIFEST.sha256'not in paths:
   need(z.getinfo('MANIFEST.sha256').file_size<=256,'checksumsidecar sizecap')
   with z.open('MANIFEST.sha256')as f:side=f.read(257)
   need(len(side)<=256,'checksumsidecar bytes');tokens=side.decode('ascii').strip().split()
   need(1<=len(tokens)<=2 and tokens[0]==manifest_sha256,'checksumsidecaridentity')
   if len(tokens)==2:need(tokens[1].lstrip('*')==root,'checksumsidecartarget')
   sidecar_verified=True;allowed_extra.add('MANIFEST.sha256')
  need(len(paths)==len(set(paths))and set(names)==set(paths)|{root}|allowed_extra,'manifest/ZIPmembership mismatch')
  if mode=='content':
   for row in rows:
    h=hashlib.sha256();count=0
    with z.open(row['path'])as f:
     for b in iter(lambda:f.read(1048576),b''):count+=len(b);need(count<=row['bytes'],'expandedbyteexcess');h.update(b)
    need(count==row['bytes']and h.hexdigest()==row['sha256'],'payloadcontentidentity')
 return dict(files=len(rows),manifest_sidecar_verified=sidecar_verified,manifest_bytes=len(raw),manifest_identity_verified=True,directory_member_sizes_verified=True,payload_content_verified=mode=='content',full_archive_SHA_recomputed=False,expanded_declared_bytes=expanded,source_provenance_paths_followed=0,numerical_execution=0,member_paths=paths)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--catalog',required=True);p.add_argument('--archive-root',required=True);p.add_argument('--mode',choices=['manifest','content'],default='manifest');p.add_argument('--report',default='evidence-check-report.json');args=p.parse_args()
 c=json.loads(Path(args.catalog).read_text());archives=c['archives'];need(0<len(archives)<=4,'boundedarchivecount');check_dependencies(archives)
 need(sum(x['manifest_bytes']for x in archives)<=CAP_MANIFEST,'totalmanifestadmissioncap');base=Path(args.archive_root).resolve();outputs=[]
 for spec in archives:
  rel=safe_relative(spec['archive']);path=(base/rel).resolve();need(path.is_relative_to(base),'archivepath escapes chosenroot');need(path.is_file()and path.stat().st_size==spec['archive_bytes'],'archive missing/sizechanged')
  r=inspect_archive(path,spec['manifest_sha256'],args.mode,spec['manifest_bytes'])
  for member in spec.get('required_entries',[]):need(safe_relative(member)in r['member_paths'],'required evidence locator missing')
  r.pop('member_paths');r.update(label=spec['label'],archive_relative_path=rel,required_evidence_entries_present=len(spec.get('required_entries',[])),historical_verified_archive_sha256=spec['historical_verified_archive_sha256']);outputs.append(r)
 result=dict(status='OFFLINE_RELATIVE_EVIDENCE_CHECK_PASS',mode=args.mode,archives=outputs,payload_content_verified_now=args.mode=='content',all_manuscript_materials_publicly_accessible=False,full_numeric_reproduction_performed=False,source_paths_followed=0,scientific_execution=0,access_status=c['access_status'])
 with Path(args.report).open('x')as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
 print(json.dumps(dict(status=result['status'],mode=args.mode,archives=len(outputs),scientific_execution=0)))
if __name__=='__main__':main()
