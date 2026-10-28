"""Pure full-native-tree covariance conversion; no parser, IO, prune or fitting."""
import math,hashlib,json
from copy import deepcopy

def need(ok,why):
 if not ok:raise ValueError(why)
def ref(r):
 need(type(r) is dict and type(r.get('path')) is str and r['path'] and type(r.get('sha256')) is str and len(r['sha256'])==64 and all(x in '0123456789abcdef' for x in r['sha256']),'explicit accepted provenance ref')
def names(values,label):
 need(type(values) is list and all(type(x) is str and x for x in values) and len(set(values))==len(values),label+' exact unique names');return list(values)
def branch(v):
 need(type(v) in (int,float) and math.isfinite(v) and v>=0,'finite nonnegative original branch length');return float(v)
def working_covariances(native_root,full54_tip_order,S_order,provenance):
 """Same full tree, original child order, explicit Root S; no mutations.
 Root provides parsed native tree and accepted identity/feasibility externally.
 The root.dist annotation is audited but contributes no incoming root edge.
 """
 full=names(full54_tip_order,'full54 tips');S=names(S_order,'S order');need(len(full)==54 and 0<len(S)<=54 and set(S)<=set(full),'exact54 full tips and explicit nonempty S subset')
 p=deepcopy(provenance)
 for k in ('source_tree_ref','rooting_acceptance_ref','root_ambiguity_ref','parser_source_ref'):ref(p[k])
 need(p['rooting_unique_claimed'] is False and type(p['candidate_root_count']) is int and p['candidate_root_count']>0 and type(p['branch_units']) is str and p['branch_units'],'explicit primary-conditional root/branch-unit provenance')
 need(native_root.up is None,'supplied tree is original full root')
 nodes=[];paths=[];depths=[];topdepths=[];objects=set();tips={};stack=[(native_root,None,0)]
 while stack:
  node,parent,level=stack.pop();need(level<=256 and len(nodes)<512,'bounded full-tree node/depth cap HOLD, no pruning')
  need(id(node) not in objects,'cycle/shared node or duplicate edge');objects.add(id(node));i=len(nodes)
  need(type(node.name) is str,'original native node name');edge=branch(node.dist);children=list(node.children)
  need((node.up is None if parent is None else node.up is nodes[parent]['_object']),'native full-tree parent/child closure')
  path=[i] if parent is None else paths[parent]+[i];depth=0. if parent is None else depths[parent]+edge;top=0 if parent is None else topdepths[parent]+1
  need(math.isfinite(depth),'finite accumulated original root path')
  nodes.append({'index':i,'parent':parent,'name':node.name,'original_dist':edge,'contributes_edge':parent is not None,'child_count':len(children),'_object':node});paths.append(path);depths.append(depth);topdepths.append(top)
  if not children:
   need(node.name and node.name not in tips,'duplicate/empty full leaf name');need(depth>0,'every full-tip original root path must be positive');tips[node.name]=i
  for child in reversed(children):stack.append((child,i,level+1))
 need(len(tips)==54 and set(tips)==set(full),'full-tree exact54 tip identity, no missing/extra/replaced names')
 def shared(i,j):
  ancestor=0
  for a,b in zip(paths[i],paths[j]):
   if a!=b:break
   ancestor=a
  return ancestor
 B=[];Bt=[];M=[]
 for a in S:
  row=[];rt=[];rm=[]
  for b in S:
   m=shared(tips[a],tips[b]);row.append(depths[m]);rt.append(float(topdepths[m]));rm.append(m)
  B.append(row);Bt.append(rt);M.append(rm)
 def normalize(matrix):
  result=[]
  for i,row in enumerate(matrix):
   di=matrix[i][i];need(di>0 and math.isfinite(di),'positive finite diagonal root path')
   out=[]
   for j,x in enumerate(row):
    dj=matrix[j][j];need(dj>0 and math.isfinite(dj),'positive finite diagonal root path')
    # Algebraic diagonal is exactly1; offdiagonal sequential division avoids
    # overflowing Bii*Bjj. This is formula evaluation, never PSD/epsilon repair.
    value=1. if i==j else x/math.sqrt(di)/math.sqrt(dj);need(math.isfinite(value),'finite normalized covariance');out.append(value)
   result.append(out)
  return result
 C=normalize(B);Ct=normalize(Bt)
 table=[{k:v for k,v in n.items() if k!='_object'} for n in nodes]
 matrix={'S_order':S,'B':B,'C':C,'B_top':Bt,'C_top':Ct}
 raw=json.dumps(matrix,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False).encode()
 return {'schema':'FULL_NATIVE54_ROOT_PATH_WORKING_COVARIANCE_V1','status':'DERIVED_UNACCEPTED_NO_MODEL_RELEASE','full54_tip_order':full,**matrix,'matrix_payload_sha256':hashlib.sha256(raw).hexdigest(),'full_node_table':table,'full_tip_paths':[{'tip':name,'node_index':tips[name],'path_node_indices':paths[tips[name]],'original_root_path':depths[tips[name]],'full_nonroot_edge_count':topdepths[tips[name]]} for name in full],'S_MRCA_node_indices':M,'source_root_name':native_root.name,'ignored_root_dist_no_incoming_edge':branch(native_root.dist),'provenance':p,'rooting_unique_claimed':False,'primary_conditional':True,'same_full_topology_no_prune_or_collapse':True,'branch_units_are_divergence_time':False,'PSD_or_epsilon_repair':False,'S_or_model_eligibility_assigned':False}
