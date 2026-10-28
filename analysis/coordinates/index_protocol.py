"""Honest schema/resource/digest protocol only; no source or window inference."""
import hashlib,json,os,sqlite3,stat,time
from pathlib import Path
SCHEMA='LOCAL_ALL_ORIGINAL_STRUCTURAL_INDEX_V1'
BIN=65536
DDL='''
CREATE TABLE source(id INTEGER PRIMARY KEY CHECK(id=1),state TEXT NOT NULL,manifest_hash TEXT NOT NULL,receipt_json BLOB);
CREATE TABLE feature_types(kind BLOB PRIMARY KEY,actual INTEGER,source_count INTEGER,emitted_count INTEGER) WITHOUT ROWID;
CREATE TABLE features(id INTEGER PRIMARY KEY,physical_line INTEGER NOT NULL,seqid BLOB NOT NULL,kind BLOB NOT NULL,raw_id BLOB NOT NULL,start0 INTEGER NOT NULL,end0 INTEGER NOT NULL,strand BLOB NOT NULL,raw_line BLOB NOT NULL);
CREATE INDEX feature_rawid ON features(raw_id,id);
CREATE TABLE bins(seqid BLOB,bin INTEGER,feature_id INTEGER,PRIMARY KEY(seqid,bin,feature_id)) WITHOUT ROWID;
CREATE INDEX bins_by_feature ON bins(feature_id,bin);
CREATE TABLE targets(ordinal INTEGER UNIQUE,target TEXT PRIMARY KEY,input_json BLOB NOT NULL);
CREATE TABLE windows(window_index INTEGER PRIMARY KEY,target TEXT,radius INTEGER,seqid_raw BLOB,requested_start0 INTEGER,requested_end0 INTEGER,RNA_start1 INTEGER,RNA_end1 INTEGER,full_original_gene_start1 INTEGER,full_original_gene_end1 INTEGER,source_gene_id_raw BLOB,strand_raw BLOB,status TEXT,gene_count_including_target INTEGER,inner_half_gene_count INTEGER,full_original_gene_length INTEGER,exon_union_bp INTEGER,union_segments INTEGER,target_gene_records INTEGER,logical_feature_window_pairs INTEGER,logical_event_keys INTEGER,primary_callable_assigned INTEGER CHECK(primary_callable_assigned=0));
CREATE TABLE exons(window_index INTEGER,segment_index INTEGER,start0 INTEGER,end0 INTEGER,PRIMARY KEY(window_index,segment_index)) WITHOUT ROWID;
CREATE TABLE typed_na(ordinal INTEGER PRIMARY KEY,input_json BLOB NOT NULL);
CREATE TABLE work_candidates(id INTEGER PRIMARY KEY);
CREATE TABLE work_ids(raw_id BLOB PRIMARY KEY,inner INTEGER) WITHOUT ROWID;
CREATE TABLE work_events(position INTEGER PRIMARY KEY,delta INTEGER);
'''

def need(ok,msg):
 if not ok:raise ValueError(msg)
def encoded(x):return json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()
def manifest_hash(m):return hashlib.sha256(encoded(m)).hexdigest()
def rowframe(fid,physical,raw):return str(fid).encode()+b'\t'+str(physical).encode()+b'\t'+str(len(raw)).encode()+b':'+raw

def limits_check(l,m):
 keys={'expected_input_bytes','expected_input_sha256','max_line_bytes','main_bytes','family_bytes','temp_bytes','owned_bytes','max_seconds','transaction_rows','max_retained_rows','max_bins_per_feature','max_bin_links','max_query_candidates','max_query_index_hits','max_feature_window_pairs','max_merged_intervals'}
 need(set(l)==keys,'exact new disk-index limits')
 need(all(type(l[k]) is int and l[k]>0 for k in keys-{'expected_input_sha256'}),'positive integer disk limits')
 need(l['expected_input_bytes']==m['expected_projection']['bytes'] and len(l['expected_input_sha256'])==64 and all(x in '0123456789abcdef' for x in l['expected_input_sha256']),'Root exact bytes/SHA; no hidden512MiB input ceiling')
 need(l['main_bytes']<=256*2**20 and l['family_bytes']<=512*2**20 and l['temp_bytes']<=128*2**20 and l['owned_bytes']<=768*2**20,'disk candidate ceilings')
 need(l['max_seconds']<=180 and l['max_line_bytes']<=2**20 and l['transaction_rows']<=1024 and l['max_bins_per_feature']<=4096 and l['max_query_candidates']<=500000,'bounded time/line/transaction/fanout/cursors')
 need(l['main_bytes']>=65536,'minimum SQLite schema space')

def file_hash(path,cap):
 p=Path(path);before=p.stat();need(stat.S_ISREG(before.st_mode) and before.st_size<=cap,'owned derived file cap')
 h=hashlib.sha256();n=0
 with p.open('rb') as f:
  while True:
   b=f.read(65536)
   if not b:break
   n+=len(b);need(n<=cap,'derived file grew');h.update(b)
 after=p.stat();need((before.st_ino,before.st_size,before.st_mtime_ns,before.st_ctime_ns)==(after.st_ino,after.st_size,after.st_mtime_ns,after.st_ctime_ns) and n==after.st_size,'closed derived file changed')
 return {'bytes':n,'sha256':h.hexdigest()}

class Guard:
 def __init__(self,db,temp,l):self.db=Path(db);self.temp=Path(temp);self.l=l;self.started=time.monotonic();self.peak={'main':0,'family':0,'temp':0};self.fault=None
 def check(self):
  need(time.monotonic()-self.started<=self.l['max_seconds'],'one cumulative180s deadline')
  main=self.db.stat().st_size if self.db.exists() else 0
  family=0
  for name in [str(self.db),str(self.db)+'-journal',str(self.db)+'-wal',str(self.db)+'-shm']:
   p=Path(name)
   try:
    need(not p.is_symlink(),'unsafe DB family symlink');family+=p.stat().st_size
   except FileNotFoundError:pass
  tmp=0
  for base,ds,fs in os.walk(self.temp,followlinks=False):
   for n in ds:need(not (Path(base)/n).is_symlink(),'unsafe temp directory symlink')
   for n in fs:
    p=Path(base)/n;need(not p.is_symlink(),'unsafe temp file symlink')
    try:tmp+=p.stat().st_size
    except FileNotFoundError:pass
  for k,v in [('main',main),('family',family),('temp',tmp)]:self.peak[k]=max(self.peak[k],v)
  need(main<=self.l['main_bytes'] and family<=self.l['family_bytes'] and tmp<=self.l['temp_bytes'] and family+tmp<=self.l['owned_bytes'],'owned DB/family/temp cap')
 def progress(self):
  try:self.check();return 0
  except BaseException as e:self.fault=e;return 1


def configure(c,temp,l,guard,*,main=True,cache_kib=8192):
 need(Path(temp).is_absolute() and Path(temp).is_dir() and not Path(temp).is_symlink(),'Root-owned explicit temp dir')
 need(type(cache_kib) is int and 0<cache_kib<=8192,'bounded SQLite cache');c.execute('PRAGMA page_size=4096');c.execute('PRAGMA cache_size=-'+str(cache_kib));c.execute('PRAGMA mmap_size=0');c.execute('PRAGMA temp_store=FILE')
 q=str(Path(temp).resolve()).replace("'","''");c.execute("PRAGMA temp_store_directory='"+q+"'")
 got=c.execute('PRAGMA temp_store_directory').fetchone();need(got and got[0]==str(Path(temp).resolve()),'explicit SQLite temp directory required; no fallback')
 c.execute('PRAGMA journal_mode=DELETE');c.execute('PRAGMA synchronous=FULL')
 c.execute('PRAGMA max_page_count='+str((l['main_bytes'] if main else l['temp_bytes'])//4096))
 c.set_progress_handler(guard.progress,1000)

def no_sort(c,query,args=()):
 for row in c.execute('EXPLAIN QUERY PLAN '+query,args):need('TEMP B-TREE' not in row[3].upper(),'unaccounted SQLite external sort forbidden')

def readonly(path,cache_kib=8192):
 c=sqlite3.connect(Path(path).resolve().as_uri()+'?mode=ro&immutable=1',uri=True);c.execute('PRAGMA query_only=ON');need(type(cache_kib) is int and 0<cache_kib<=8192,'bounded readonly cache');c.execute('PRAGMA cache_size=-'+str(cache_kib));c.execute('PRAGMA mmap_size=0');return c
