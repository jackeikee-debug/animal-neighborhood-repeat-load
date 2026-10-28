"""One new complete G or S annotation-union projection, no raw persistence.

Raw reference annotation is permitted input, not held-out CNE outcomes. All
forbidden source portions are removed before class/union use; rows entirely
forbidden are not class-parsed. Such holes remain explicitly unknown, and every
released inverse target footprint must avoid them. No protected burden output.
"""
from pathlib import Path
from array import array
from bisect import bisect_right
import gzip,hashlib,json,re,shlex,sys,urllib.request
import numpy as np
from physical_shift import AnnotationUnion,prepare_callable_fragments,query_counts
R3=Path(__import__('os').environ['ANRL_R3_ROOT'])
P=Path(__import__('os').environ['ANRL_PROJECT_ROOT'])
cfg=json.loads(Path(sys.argv[1]).read_text());sp=cfg['species'];assert sp in ['G','S']
geo=json.loads(gzip.decompress((R3/'geometry-design-preparation-syntaxfix-once/geometry-and-mask-design.json.gz').read_bytes()))
s=geo['species'][sp];ids=geo['ids'];N=len(ids);assert N==1896
assert array('I').itemsize==4 and max(s['lengths'].values())<2**32, 'uint32_source_buffer_contract'
source=json.loads((R3/'G-S-known-source-metadata-first/accepted-source-identities.json').read_text())['sources'][sp]
env={}
for line in(P/'galaxy/galaxy.env').read_text().splitlines():
    fields=shlex.split(line.removeprefix('export '),comments=True)
    if fields and '=' in fields[0]:k,v=fields[0].split('=',1);env[k]=v
class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self,*a,**k):raise RuntimeError('source_redirect_forbidden')
opener=urllib.request.build_opener(NoRedirect)
blocked={ch:(z,[b for a,b in z])for ch,z in s['forbidden'].items()}
def outside_forbidden(ch,a,b):
    z,ends=blocked.get(ch,([],[]));i=bisect_right(ends,a);out=[];cursor=a
    while i<len(z)and z[i][0]<b:
        left,right=z[i]
        if cursor<left:out.append((cursor,min(left,b)))
        cursor=max(cursor,right)
        if cursor>=b:break
        i+=1
    if cursor<b:out.append((cursor,b))
    return out
buffers={ch:array('I')for ch in s['chromosomes']};seen_known=set();NR=0;bytecount=0;digest=hashlib.sha256();included=0
number=re.compile(r'^[0-9]+(?:\.[0-9]+)?$');unsigned=re.compile(r'^\(?[0-9,]+\)?$');leftnum=re.compile(r'^\(?-?[0-9,]+\)?$');classes={'DNA','LINE','LTR','SINE','RC','Retroposon','PLE','SVA','Unknown'}
header='SW score\t% div.\t% del.\t% ins.\tquery sequence\tpos in  query: begin\tend\t(left)\trepeat\tclass/family\tpos in repeat: begin\tend\t(left)\tID'
def consume(raw):
    global NR,included
    NR+=1;assert len(raw)<16384
    text=raw.decode().rstrip('\r\n')
    if NR==1:assert text==header;return
    if NR==2:assert text=='';return
    f=text.split('\t');star=f[-1]=='*';n=len(f)-int(star)
    assert n in (14,15)and f[4]in s['lengths'] and f[5].isdigit()and f[6].isdigit()
    ch=f[4];a=int(f[5])-1;b=int(f[6]);assert 0<=a<b<=s['lengths'][ch]
    if str(NR)in s['known_exceptions']:
        assert text==s['known_exceptions'][str(NR)];seen_known.add(str(NR));return
    if ch not in buffers:return
    pieces=outside_forbidden(ch,a,b)
    if not pieces:return
    assert all(number.fullmatch(x)for x in f[:4])and re.fullmatch(r'\([0-9,]+\)',f[7])
    if f[8]=='C':assert n==15;rep,cl,rl,re_,rb,rid=f[9:15]
    else:assert n==14;rep,cl,rb,re_,rl,rid=f[8:14]
    assert rep and cl and rid.isdigit()and unsigned.fullmatch(rb)and unsigned.fullmatch(re_)and leftnum.fullmatch(rl)
    if cl.split('/')[0]not in classes:return
    included+=1
    for l,r in pieces:buffers[ch].extend((l,r))
url=env['GALAXY_URL'].rstrip('/')+source['download_route']
Path('source-body-intent-once.json').write_text(json.dumps(dict(species=sp,source=source,expected_bytes=source['file_size'],expected_rows=s['expected_source_rows'],new_source_projection=True,RepeatMasker_execution=0,raw_persistence=False,automatic_retry=False))+'\n')
with opener.open(urllib.request.Request(url,headers={'x-api-key':env['GALAXY_API_KEY'],'Accept-Encoding':'identity'}),timeout=45)as r:
    assert r.status==200 and r.headers.get('Content-Encoding')in(None,'identity')
    if r.headers.get('Content-Length')is not None:assert int(r.headers['Content-Length'])==source['file_size']
    tail=b''
    while True:
        chunk=r.read(2**20)
        if not chunk:break
        bytecount+=len(chunk);assert bytecount<=source['file_size'];digest.update(chunk)
        z=(tail+chunk).split(b'\n');tail=z.pop();assert len(tail)<16384, 'unterminated_source_line_cap'
        for line in z:consume(line)
    if tail:consume(tail)
assert bytecount==source['file_size']and NR==s['expected_source_rows']+2 and seen_known==set(s['known_exceptions'])
payload={};meta=[];zero=np.zeros(N,dtype=np.int64);positions={g:i for i,g in enumerate(ids)}
for j,ch in enumerate(s['chromosomes']):
    raw=np.frombuffer(buffers[ch],dtype=np.uint32).reshape(-1,2);assert len(raw)>0
    order=np.argsort(raw[:,0],kind='stable');ordered=raw[order];merged=array('I')
    a,b=map(int,ordered[0])
    for l,r in ordered[1:]:
        l,r=int(l),int(r)
        if l<=b:b=max(b,r)
        else:merged.extend((a,b));a,b=l,r
    merged.extend((a,b));z=np.frombuffer(merged,dtype=np.uint32).reshape(-1,2)
    starts=z[:,0].astype(np.int64);ends=z[:,1].astype(np.int64);cumulative=np.r_[np.int64(0),np.cumsum(ends-starts,dtype=np.int64)]
    L=s['lengths'][ch];A=AnnotationUnion(starts,ends,cumulative,L)
    targets={t['orthogroup']:t for t in s['targets']if t['seqid']==ch}
    fr=prepare_callable_fragments({g:targets[g]['callable']if g in targets else []for g in ids},L,target_ids=ids)
    counts=query_counts(A,fr,0,L,circular=True);zero+=counts
    payload[f'c{j}_starts']=starts;payload[f'c{j}_ends']=ends;payload[f'c{j}_cumulative']=cumulative
    meta.append(dict(chromosome=ch,L=L,key=f'c{j}',interval_count=len(starts),retained_union_bp=int(cumulative[-1]),forbidden_holes_explicit=True))
    del buffers[ch],raw,order,ordered,merged,z,A
expected=np.array([t['A']for t in s['targets']],dtype=np.int64)
assert np.array_equal(zero,expected),('zero_phase_numerator_not_identical',sp,int(np.sum(zero!=expected)))
np.savez_compressed('accepted-A-unions.npz',**payload)
Path('track-metadata.json').write_text(json.dumps(dict(species=sp,chromosomes=meta,forbidden=s['forbidden'],source=source,source_body_bytes=bytecount,source_body_sha256=digest.hexdigest(),source_rows=NR-2,known_exceptions=sorted(seen_known,key=int),original_zero_phase_numerators_exact=1896,old_models_recomputed=0,raw_source_persistence=0,annotation_producer_recomputed=0,protected_outcomes_opened=0),indent=2)+'\n')
Path('return.json').write_text(json.dumps(dict(status='NEW_COMPLETE_NEEDED_CHROMOSOME_A_TRACK_AND_ALL1896_ZERO_PHASE_COUNTS_ACCEPTED',species=sp,source_body_bytes=bytecount,source_rows=NR-2,intervals=sum(x['interval_count']for x in meta),all_zero_phase_numerators_exact=True,real_null_draws=0,old_fits=0,raw_persistence=0),indent=2)+'\n')
