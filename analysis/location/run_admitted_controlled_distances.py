"""Root's one finite first evaluation of unique admitted controlled offsets."""
from pathlib import Path
import csv,gzip,hashlib,json,sys
import numpy as np
from location_binding import validate_freeze,validate_draw_rows,phase_matches,validate_pair
Q4=Path(__import__('os').environ['ANRL_R4_ROOT']);R=Q4.parent
sys.path.insert(0,str(R/'code'))
from physical_shift_root_numeric_gates import AnnotationUnion,prepare_callable_fragments,query_counts,fixed_qr,response_residuals,pearson
contract=json.loads(Path(sys.argv[1]).read_text())
contract_pins={str(Path(k).resolve()):v for k,v in contract['pins'].items()}
def checkpin(path):
 path=Path(path);assert hashlib.sha256(path.read_bytes()).hexdigest()==contract_pins[str(path.resolve())],('Rootbinding',path)
freeze_path=Q4/'geometry-and-saved1998-first/fixed-distance-admission-freeze.json';checkpin(freeze_path)
freeze=json.loads(freeze_path.read_text());validate_freeze(freeze,contract['expected_unique_keys']);cases=freeze['unique_admitted_cases']
assert freeze['status']=='GEOMETRY_ONLY_FIXED72_SCENARIOS_FROZEN_BEFORE_NEW_RESPONSES'and len(freeze['cases'])==72
gfile=R/'geometry-design-preparation-syntaxfix-once/geometry-and-mask-design.json.gz';raw=gfile.read_bytes();geo=json.loads(gzip.decompress(raw));bind=json.loads((R/'geometry-design-preparation-syntaxfix-once/Root-paired-geometry-X-row-order-binding.json').read_text())
assert hashlib.sha256(raw).hexdigest()==bind['geometry_sha256']and geo['ids']==bind['oldgroup_ids'];ids=geo['ids'];N=len(ids);assert N==1896
zfile=R/'geometry-design-preparation-syntaxfix-once/original-fixed-design-and-saved-residuals.npz';assert hashlib.sha256(zfile.read_bytes()).hexdigest()==bind['npz_sha256'];z=np.load(zfile,allow_pickle=False)
P=z['P_fixed_residual'];assert P.shape==(N,)and np.isfinite(P).all();basis={sp:fixed_qr(z[sp+'_X'])for sp in ['G','S']}
phase=R/'all1998-offsets-geometry-freeze-first';pm=json.loads((phase/'return.json').read_text());plans=json.loads((phase/'allowed-phase-plans.json').read_text())
for f in [phase/'return.json',phase/'allowed-phase-plans.json',Q4/'Root-geometry-and-saved-statistic-identity-acceptance.json']:checkpin(f)
accepted=json.loads((R/'physical-output-complete-QA-first/accepted-two-mode-physical-backgrounds.json').read_text())
for mode in ['circ','linear']:
 f=phase/(mode+'-all999-frozen-offsets.npy');checkpin(f);pin=next(x for x in pm['outputs']if x['mode']==mode);assert hashlib.sha256(f.read_bytes()).hexdigest()==pin['sha256']
 for name in ['all999-conditional-background-statistics.tsv','all999-two-species-frozen-window-counts.npz']:
  f=R/(mode+'-all999-first')/name;checkpin(f);pin=next(x for x in accepted['outputs']if Path(x['path']).parent.name==mode+'-all999-first'and Path(x['path']).name==name);assert f.stat().st_size==pin['bytes']and hashlib.sha256(f.read_bytes()).hexdigest()==pin['sha256']
old_offsets={mode:np.load(phase/(mode+'-all999-frozen-offsets.npy'),allow_pickle=False)for mode in ['circ','linear']};old_rows={}
for mode in ['circ','linear']:
 with(R/(mode+'-all999-first')/'all999-conditional-background-statistics.tsv').open()as f:old_rows[mode]=list(csv.DictReader(f,delimiter='\t'))
 assert old_offsets[mode].shape==(999,45)and old_offsets[mode].dtype==np.dtype('int64');validate_draw_rows(old_rows[mode])
tracks={};denom={};columns=[]
for sp in ['G','S']:
 s=geo['species'][sp];stage=R/(sp+'-original-annotation-projection-first');m=json.loads((stage/'track-metadata.json').read_text());proof=json.loads((stage/'Root-track-acceptance.json').read_text())
 assert m['species']==sp and m['original_zero_phase_numerators_exact']==1896
 for pin in proof['outputs']:
  f=stage/Path(pin['path']).name
  assert f.stat().st_size==pin['bytes']and hashlib.sha256(f.read_bytes()).hexdigest()==pin['sha256']
 payload=np.load(stage/'accepted-A-unions.npz',allow_pickle=False);denom[sp]=np.array([t['D']for t in s['targets']],dtype=np.int64);assert[t['orthogroup']for t in s['targets']]==ids
 parts=[];seen=set();total=np.zeros(N,dtype=np.int64)
 for cm in m['chromosomes']:
  ch=cm['chromosome'];L=cm['L'];key=cm['key'];assert L==s['lengths'][ch]
  A=AnnotationUnion(payload[key+'_starts'],payload[key+'_ends'],payload[key+'_cumulative'],L)
  assert np.all(A.starts<A.ends)and np.all(A.ends[:-1]<A.starts[1:])and A.cumulative[0]==0 and np.array_equal(np.diff(A.cumulative),A.ends-A.starts)
  targets={t['orthogroup']:t for t in s['targets']if t['seqid']==ch};assert not seen.intersection(targets);seen.update(targets)
  fr=prepare_callable_fragments({g:targets[g]['callable']if g in targets else[]for g in ids},L,target_ids=ids,expected_D=[targets[g]['D']if g in targets else 0 for g in ids]);total+=fr.denominators
  parts.append((A,fr,L,ch));columns.append([sp,ch])
 assert seen==set(ids)and np.array_equal(total,denom[sp]);tracks[sp]=parts
assert len(columns)==45 and all(x['columns']==columns for x in pm['outputs'])
lengths=np.array([geo['species'][sp]['lengths'][ch]for sp,ch in columns],dtype=np.int64)
results=[];saved_counts=[];projected=0;old_reuse=0;progress=Path('unique-scenario-progress.jsonl').open('x')
for case in cases:
 key=case['canonical_key'];mode=key.split(':')[0];assert mode in ['circ','linear']and case['admissible']and (case['G_delta']or case['S_delta'])
 vector=np.array([case[sp+'_delta']%geo['species'][sp]['lengths'][ch]if mode=='circ'else case[sp+'_delta']for sp,ch in columns],dtype=np.int64)
 matches=[]
 for oldmode in ['circ','linear']:
  # Both accepted operators query the same modulo phase when a no-wrap
  # realization exists; compare all prior vectors across mode labels.
  rows=phase_matches(old_offsets[oldmode],vector,lengths)
  if len(rows):matches.append((oldmode,int(rows[0])))
 progress.write(json.dumps(dict(canonical_key=key,state='STARTED',old_saved_match=matches))+'\n');progress.flush()
 if matches:
  oldmode,index=matches[0];row=old_rows[oldmode][index];rec={k:float(row[k])for k in ['r_GP','r_GS','r_PS','T','G_residual_norm','S_residual_norm']}
  with np.load(R/(oldmode+'-all999-first')/'all999-two-species-frozen-window-counts.npz',allow_pickle=False)as old:pair=np.array(old['counts'][index],copy=True)
  validate_pair(pair,np.stack([denom['G'],denom['S']]))
  old_reuse+=1;origin='EXACT_OLD_SAVED_OFFSET_VECTOR_REUSE'
 else:
  residuals={};pair=np.empty((2,N),dtype=np.uint32)
  for spi,sp in enumerate(['G','S']):
   counts=np.zeros(N,dtype=np.int64);delta=case[sp+'_delta']
   for A,fr,L,ch in tracks[sp]:
    value=delta%L if mode=='circ'else delta
    assert any(a<=value<b for a,b in plans[mode][sp][ch]['allowed']),('unreleased point',key,sp,ch)
    counts+=query_counts(A,fr,value,L,circular=mode=='circ')
   assert np.all(counts>=0)and np.all(counts<=denom[sp]);pair[spi]=counts
   residuals[sp]=response_residuals(counts,denom[sp],basis[sp]);projected+=1
  GP=pearson(residuals['G'],P);GS=pearson(residuals['G'],residuals['S']);PS=pearson(P,residuals['S'])
  rec=dict(r_GP=GP,r_GS=GS,r_PS=PS,T=float((GP+GS+PS)/3),G_residual_norm=float(np.linalg.norm(residuals['G'])),S_residual_norm=float(np.linalg.norm(residuals['S'])));origin='NEW_CONTROLLED_OFFSET_FIRST_EVALUATION'
 validate_pair(pair,np.stack([denom['G'],denom['S']]));assert np.isfinite(list(rec.values())).all()
 rec.update(canonical_key=key,distance_bp=case['distance_bp'],G_delta=case['G_delta'],S_delta=case['S_delta'],canonical_mode=mode,origin=origin,old_origin_mode=matches[0][0]if matches else'',old_origin_draw=matches[0][1]+1 if matches else'',count_index=len(saved_counts));results.append(rec);saved_counts.append(pair)
 progress.write(json.dumps(dict(canonical_key=key,state='RETURNED',T=rec['T']))+'\n');progress.flush()
progress.close();assert len(results)==len(cases)
with Path('unique-controlled-offset-results.tsv').open('x')as f:
 w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(results)
counts=np.stack(saved_counts);np.savez_compressed('unique-controlled-offset-counts.npz',counts=counts)
Path('return.json').write_text(json.dumps(dict(status='ALL_ADMITTED_UNIQUE_CONTROLLED_OFFSETS_FIRST_EVALUATION_COMPLETE',logical_scenarios=72,logical_admitted=sum(c['admissible']for c in freeze['cases']),unique_results=len(results),new_rank_QR_projections=projected,old_saved_vector_reuses=old_reuse,old_models_refitted=0,new_source_reads=0,old_nulls_reexecuted=0,new_pvalues=0,per_species_members=N,counts_shape=list(counts.shape),count_array_sha256=hashlib.sha256(counts.tobytes()).hexdigest(),undefined_points=0,omissions=0,scope='Only admitted near-zero genomic offsets; no longrange specificity claim'),indent=2)+'\n')
