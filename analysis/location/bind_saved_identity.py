"""New compact identity gate closes I1 without rerunning geometry/statistics."""
from pathlib import Path
import csv,hashlib,json
import numpy as np
from location_binding import validate_freeze,validate_draw_rows
Q=Path(__import__('os').environ['ANRL_R4_ROOT']);R=Q.parent
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def desc(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p))
freeze_path=Q/'geometry-and-saved1998-first/fixed-distance-admission-freeze.json';freeze=json.loads(freeze_path.read_text());keys=sorted(x['canonical_key']for x in freeze['unique_admitted_cases']);validate_freeze(freeze,keys)
accepted=json.loads((R/'physical-output-complete-QA-first/accepted-two-mode-physical-backgrounds.json').read_text());pm=json.loads((R/'all1998-offsets-geometry-freeze-first/return.json').read_text());bindings=[]
with(Q/'geometry-and-saved1998-first/saved1998-distance-and-footprint-T.tsv').open()as f:joined=list(csv.DictReader(f,delimiter='\t'))
for mode in ['circ','linear']:
 stage=R/(mode+'-all999-first')
 for name in ['all999-conditional-background-statistics.tsv','all999-two-species-frozen-window-counts.npz']:
  p=stage/name;record=next(x for x in accepted['outputs']if Path(x['path']).parent.name==stage.name and Path(x['path']).name==name);assert p.stat().st_size==record['bytes']and sha(p)==record['sha256'];bindings.append(desc(p))
 with(stage/'all999-conditional-background-statistics.tsv').open()as f:rows=list(csv.DictReader(f,delimiter='\t'))
 validate_draw_rows(rows);j=[x for x in joined if x['mode']==mode];validate_draw_rows(j)
 assert all(float(a['T'])==float(b['T_reused'])for a,b in zip(rows,j))
 p=R/'all1998-offsets-geometry-freeze-first'/(mode+'-all999-frozen-offsets.npy');m=next(x for x in pm['outputs']if x['mode']==mode);assert sha(p)==m['sha256'];a=np.load(p,allow_pickle=False);assert a.shape==(999,45)and a.dtype==np.dtype('int64');bindings.append(desc(p))
report=dict(status='ROOT_ACCEPTED_I1_FREEZE_AXIS_AND_SAVED1998_STATISTIC_IDENTITIES_NO_REPLAY',freeze=desc(freeze_path),expected_unique_keys=keys,logical_cases=72,admitted_aliases=16,unique=8,saved_original_output_bindings=bindings,saved_join=desc(Q/'geometry-and-saved1998-first/saved1998-distance-and-footprint-T.tsv'),all1998_T_equal_accepted_saved_rows=True,numerical_replays=0,NEW_review_binding_fixtures=8)
Path('Root-geometry-and-saved-statistic-identity-acceptance.json').write_text(json.dumps(report,indent=2)+'\n')
