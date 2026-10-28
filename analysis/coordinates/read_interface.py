"""Read-only bounded window/exon interface, no aliases/Genome/eligibility."""
import hashlib,json
from index_protocol import SCHEMA,need,file_hash,readonly,encoded

class View:
 def __init__(self,path,proof):
  need(proof['schema']==SCHEMA and proof['status']=='DERIVED_INDEX_VALIDATED_SOURCE_ONLY' and proof['primary_callable_or_model_eligible'] is False,'independently validated exact index required')
  need(file_hash(path,256*2**20)==proof['db_file'],'pinned closed derived DB')
  self.c=readonly(path)
  try:
   r=self.c.execute('SELECT state,manifest_hash,receipt_json FROM source').fetchone()
   need(r and r[0]=='CLOSED' and r[1]==proof['manifest_sha256'] and hashlib.sha256(r[2]).hexdigest()==proof['source_receipt_sha256'],'validated source ledger pin')
  except BaseException:self.c.close();raise
 def __enter__(self):return self
 def __exit__(self,*args):self.c.close()
 def windows(self):
  names=[r[1] for r in self.c.execute('PRAGMA table_info(windows)')]
  for row in self.c.execute('SELECT * FROM windows ORDER BY window_index'):
   r=dict(zip(names,row));r['target_key']=r['target'];yield r
 def exons(self,window_index):
  need(type(window_index) is int and window_index>0 and self.c.execute('SELECT 1 FROM windows WHERE window_index=?',(window_index,)).fetchone() is not None,'exact existing window index')
  for row in self.c.execute('SELECT start0,end0 FROM exons WHERE window_index=? ORDER BY segment_index',(window_index,)):yield row
 def typed_na(self):
  for (raw,) in self.c.execute('SELECT input_json FROM typed_na ORDER BY ordinal'):yield json.loads(raw)
 def source_receipt(self):return json.loads(self.c.execute('SELECT receipt_json FROM source').fetchone()[0])

def open_view(db_path,validated_receipt):return View(db_path,validated_receipt)
