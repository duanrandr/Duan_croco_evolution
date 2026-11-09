"""Recompute the original three-trait mean-profile analysis for 25 treatments."""
from pathlib import Path
import argparse, csv, hashlib, json
import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform

# Step 1: Set up treatment order and analysis definitions.
ROOT=Path(__file__).resolve().parents[2]
TRAITS=['Growth_rates','Nfix_PON','Cfix_POC']
GROUPS=['T0_28HFe']
for temp, transfers in [(23,[28,32]),(28,[23,32]),(32,[23,28])]:
    for prefix in [f'IA_{temp}',str(temp)]+[f'{temp}-{end}'for end in transfers]:
        GROUPS.extend(prefix+fe for fe in ['HFe','LFe'])
TARGET={'23LFe','23-28LFe','23-32LFe'}

def save(name,rows):
    with (OUT/name).open('w',newline='')as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

# Step 2: Read canonical measurements and choose a separate output location.
parser=argparse.ArgumentParser();parser.add_argument('--raw',required=True);parser.add_argument('--summary',required=True)
parser.add_argument('--out',default=str(ROOT/'results/reproduced/multivariate'))
args=parser.parse_args();raw_path=Path(args.raw).resolve();summary_path=Path(args.summary).resolve()
OUT=Path(args.out)/'source_data'; OUT.mkdir(parents=True,exist_ok=True)
raw=list(csv.DictReader(raw_path.open()));summary=list(csv.DictReader(summary_path.open()))
assert len(raw)==75 and len(summary)==75,(len(raw),len(summary))
assert set(r['Treatment_Group']for r in raw)==set(GROUPS)
assert set(r['Treatment_Group']for r in summary)==set(GROUPS)
lookup={(r['Treatment_Group'],r['Trait']):r for r in summary}
assert len(lookup)==75
matrix=np.array([[float(lookup[g,t]['Mean']) for t in TRAITS]for g in GROUPS])
assert matrix.shape==(25,3)and np.isfinite(matrix).all()
errors=[];ns={}
for i,g in enumerate(GROUPS):
    for j,t in enumerate(TRAITS):
        vals=[float(r[t])for r in raw if r['Treatment_Group']==g and r[t]not in ('','NA','NaN','nan','None')]
        ns[g,t]=len(vals);assert len(vals)==int(float(lookup[g,t]['n']))
        assert len(vals)==(2 if(g,t)==('23-32HFe','Cfix_POC')else 3)
        errors.append(abs(np.mean(vals)-matrix[i,j]))
assert max(errors)<1e-12
# Step 3: Standardize treatment means and compute Euclidean distances.
means=matrix.mean(axis=0);sd=matrix.std(axis=0,ddof=1);z=(matrix-means)/sd
assert np.max(np.abs(z.mean(axis=0)))<1e-12
assert np.max(np.abs(z.std(axis=0,ddof=1)-1))<1e-12
metadata=[]
for i,g in enumerate(GROUPS):
    fe='HFe'if g.endswith('HFe')else 'LFe';stem=g.removesuffix(fe)
    if stem.startswith('T0_'):
        stage='T0';short='Ancestral (T0)\n(28°C HFe)';long='Ancestral (T0) | 28°C | HFe'
    elif stem.startswith('IA_'):
        stage='IA';temp=stem[3:];short=f'IA {temp}{fe}';long=f'Initial acclimation | {temp}°C | {fe}'
    elif '-'in stem:
        stage='Transfer';short=stem.replace('-','→')+fe;long=f'Transfer | {stem.replace("-","→")}°C | {fe}'
    else:
        stage='TF';short=f'TF {stem}{fe}';long=f'TF | {stem}°C | {fe}'
    metadata.append({'Treatment_Group':g,'Stage':stage,'Fe':fe,'Highlight':g in TARGET,'Heatmap_Label':short,'Dendrogram_Label':long,**{t:float(matrix[i,j])for j,t in enumerate(TRAITS)},**{t+'_n':ns[g,t]for t in TRAITS}})
save('Treatment_mean_matrix.csv',metadata)
save('Summary_All_canonical.csv',summary)
save('Raw_profiles_reference.csv',raw)
save('Standardized_matrix.csv',[{'Treatment_Group':g,**{t:float(z[i,j])for j,t in enumerate(TRAITS)}}for i,g in enumerate(GROUPS)])
save('Standardization_parameters.csv',[{'Trait':t,'Grand_mean_of_25_treatment_means':float(means[j]),'Sample_SD_of_25_treatment_means':float(sd[j]),'ddof':1,'n_treatments':25}for j,t in enumerate(TRAITS)])
dmat=squareform(pdist(z,metric='euclidean'))
save('Euclidean_distance_matrix.csv',[{'Treatment_Group':g,**{other:float(dmat[i,j])for j,other in enumerate(GROUPS)}}for i,g in enumerate(GROUPS)])

# Step 4: Preserve the prior UPGMA procedure and child-index orientation.
n=len(GROUPS);active={i:[i]for i in range(n)}
distance={(i,j):float(dmat[i,j])for i in range(n)for j in range(i+1,n)}
rows=[];members={i:[i]for i in range(n)}
for step in range(n-1):
    a,b=min(distance,key=lambda key:(distance[key],key));height=distance[a,b]
    na,nb=len(active[a]),len(active[b]);node=n+step
    updates={tuple(sorted((k,node))):(na*distance[tuple(sorted((a,k)))]+nb*distance[tuple(sorted((b,k)))])/(na+nb)for k in active if k not in(a,b)}
    distance={key:value for key,value in distance.items()if a not in key and b not in key};distance.update(updates)
    children=active.pop(a)+active.pop(b);active[node]=children;members[node]=children
    rows.append({'Node_ID':node,'Child_1':a,'Child_2':b,'Height':height,'Count':len(children),'All_highlight':all(GROUPS[k]in TARGET for k in children)})
save('Average_linkage.csv',rows)
order=members[2*n-2]
save('Leaf_order_bottom_to_top.csv',[{'Position':p+1,'Input_Index':i,'Treatment_Group':GROUPS[i]}for p,i in enumerate(order)])

# Step 5: Verify and save the independent SciPy comparison by clade, irrespective of child orientation.
sp=linkage(z,method='average',metric='euclidean');sm={i:[i]for i in range(n)};sh={}
for k,(a,b,h,c)in enumerate(sp):
    inds=sm[int(a)]+sm[int(b)];sm[n+k]=inds;sh[tuple(sorted(inds))]=float(h)
validation=[]
for row in rows:
    key=tuple(sorted(members[row['Node_ID']]))
    assert key in sh
    err=abs(row['Height']-sh[key]);assert err<1e-12
    validation.append({'Node_ID':row['Node_ID'],'Manual_height':row['Height'],'SciPy_height':sh[key],'Absolute_difference':err,'Members':';'.join(GROUPS[k]for k in key)})
save('SciPy_linkage_validation.csv',validation)
highlight=[{'Node_ID':r['Node_ID'],'Height':r['Height'],'Members':[GROUPS[k]for k in members[r['Node_ID']]]}for r in rows if r['All_highlight']]
report={'input_raw':raw_path.name,'input_raw_sha256':hashlib.sha256(raw_path.read_bytes()).hexdigest(),'input_summary':summary_path.name,'input_summary_sha256':hashlib.sha256(summary_path.read_bytes()).hexdigest(),'matrix_shape':[25,3],'stages':{'T0':1,'IA':6,'TF':6,'Transfer':12},'analysis_unit':'Treatment-group mean profile','traits':TRAITS,'scaling':'z scores across all 25 treatment means for each variable; sample SD (ddof=1)','clustering':'Euclidean distance, average linkage (UPGMA)','raw_to_summary_max_error':float(max(errors)),'missing_policy':'No imputation; observed means per trait and treatment. 23-32HFe Cfix_POC n=2, all other treatment-by-trait n=3.','target_groups':sorted(TARGET),'actual_pure_target_clades':highlight,'three_target_groups_form_pure_clade':any(len(r['Members'])==3 for r in highlight),'max_SciPy_height_error':max(r['Absolute_difference']for r in validation),'all_true_T0_traits_present':True,'new_true_T0_profile':dict(zip(TRAITS,map(float,matrix[0]))),'no_stage_collapsing':True}
(OUT/'matrix_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
