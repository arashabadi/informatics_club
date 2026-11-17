from pathlib import Path
import numpy as np,json,importlib.metadata
from scipy.spatial.distance import cdist
import umap
from umap.umap_ import smooth_knn_dist

root=Path(__file__).resolve().parents[1]
cells=list('ABCDEF');genes=['CD3D','MS4A1','NKG7','LST1']
counts=np.array([[4,0,1,0],[8,0,2,0],[0,6,0,0],[1,10,1,0],[0,0,12,8],[1,0,16,10]],dtype=float)
totals=np.array([1000,2000,1200,2000,1500,2000],dtype=float)
normalized=counts/totals[:,None]*1e4
y=np.log1p(normalized);means=y.mean(0);sd=y.std(0,ddof=1);z=(y-means)/sd
u,s,vt=np.linalg.svd(z,full_matrices=False)
# Deterministic sign convention: the largest-magnitude loading is positive.
for i in range(vt.shape[0]):
 if vt[i,np.argmax(abs(vt[i]))]<0:vt[i]*=-1
loadings=vt.T;scores=z@loadings;ratio=s*s/np.sum(s*s);pca=scores[:,:3]
dist=cdist(pca,pca)
model=umap.UMAP(n_neighbors=3,n_components=2,min_dist=.2,metric='euclidean',
 unique=True,random_state=42,n_jobs=1,n_epochs=500,init='random')
embedding=model.fit_transform(pca)
_,idx,inv=np.unique(pca,axis=0,return_index=True,return_inverse=True)
up=pca[idx];ud=cdist(up,up)
nearest=np.argsort(ud,axis=1,kind='stable')[:,:3]
knnd=np.take_along_axis(ud,nearest,axis=1)
sigma,rho=smooth_knn_dist(knnd.astype(np.float32),3,local_connectivity=1)
directed=np.zeros_like(ud)
for i in range(len(idx)):
 for j in nearest[i,1:]:directed[i,j]=np.exp(-max(0,ud[i,j]-rho[i])/sigma[i])
sym=directed+directed.T-directed*directed.T
np.fill_diagonal(sym,0)
assert np.allclose(sym,model.graph_.toarray(),atol=3e-5), (sym,model.graph_.toarray())
np.testing.assert_allclose(pca[0],pca[1],atol=1e-10)
np.testing.assert_allclose(embedding[0],embedding[1],atol=1e-10)
np.testing.assert_allclose(cdist(z,z),cdist(scores,scores),atol=1e-12)
np.testing.assert_allclose(z.mean(0),0,atol=1e-12)
np.testing.assert_allclose(z.std(0,ddof=1),1,atol=1e-12)
labels=['A/B' if j in [0,1] else cells[j] for j in idx]
# Full pairwise cross-entropy is for teaching the target, not the sampled fit loss.
a=float(model._a);b=float(model._b)
initial=np.random.RandomState(17).uniform(-1,1,(len(idx),2))
fin=np.array([embedding[j] for j in idx])
q=lambda yy:1/(1+a*cdist(yy,yy)**(2*b))
qi=q(initial);qf=q(fin)
def ce(p,q):
 q=np.clip(q,1e-8,1-1e-8)
 return -p*np.log(q)-(1-p)*np.log1p(-q)
pairs=[]
for i in range(len(idx)):
 for j in range(i+1,len(idx)):
  pairs.append(dict(pair=labels[i]+' / '+labels[j],i=i,j=j,p=float(sym[i,j]),
   distance_PC=float(ud[i,j]),q_initial=float(qi[i,j]),q_final=float(qf[i,j]),
   ce_initial=float(ce(sym[i,j],qi[i,j])),ce_final=float(ce(sym[i,j],qf[i,j]))))
data=dict(cells=cells,genes=genes,counts=counts.tolist(),totals=totals.tolist(),normalized=normalized.tolist(),
 log1p=y.tolist(),mean=means.tolist(),sd=sd.tolist(),z=z.tolist(),loadings=loadings.tolist(),scores=scores.tolist(),
 variance_ratio=ratio.tolist(),eigenvalues=(s*s/5).tolist(),correlation=np.corrcoef(z.T).tolist(),
 distances=dist.tolist(),unique_indices=idx.tolist(),inverse=inv.tolist(),graph_labels=labels,
 unique_distances=ud.tolist(),nearest=nearest.tolist(),knn_distances=knnd.tolist(),sigma=sigma.tolist(),rho=rho.tolist(),
 directed=directed.tolist(),graph=sym.tolist(),embedding=embedding.tolist(),a=a,b=b,pairs=pairs,
 versions={p:importlib.metadata.version(p) for p in ['numpy','scipy','scikit-learn','umap-learn']},
 parameters=dict(n_neighbors=3,n_components=2,min_dist=.2,metric='euclidean',unique=True,random_state=42,n_epochs=500,init='random',n_jobs=1))
(root/'build/numerical-example.json').write_text(json.dumps(data,indent=2))
print('variance ratio',ratio.round(6));print('scores',scores.round(4));print('labels',labels)
print('rho',rho);print('sigma',sigma);print('knn',nearest);print('graph',sym.round(4));print('embedding',embedding.round(4));print('a,b',a,b);print('pairs',pairs)
