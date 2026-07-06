"""Regenerate mos_field.png and mos_ppc.png from the NUTS posterior chains."""
import numpy as np, glob
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
rng=np.random.default_rng(0)
d=np.load("mos_arrays.npz",allow_pickle=True)
YC,YP=d['YC'],d['YP']; m,T=YC.shape
elev=np.nan_to_num((d['elev']-np.nanmean(d['elev']))/(np.nanstd(d['elev'])+1e-9))
wet=d['wetland'].astype(float); forr=d['forest'].astype(float)
W=d['W'].astype(float); deg=np.maximum(W.sum(1),1.0)
obsC=np.isfinite(YC); obsP=np.isfinite(YP)
ch=[np.load(f) for f in sorted(glob.glob("mos_chain_*.npz"))]
def cat(k): return np.concatenate([c[k] for c in ch])
z0=cat('z0'); eps=cat('eps'); r=cat('r'); kap=cat('kappa'); eta=cat('eta')
phi=cat('phi'); psi=cat('psi'); delta=cat('delta'); nu=cat('nu'); sig=cat('sigma_p')
b0=cat('b0'); be=cat('b_elev'); bw=cat('b_wet'); bf=cat('b_for')
S=z0.shape[0]
# thin to 400 draws for speed
idx=rng.choice(S,size=min(400,S),replace=False)
def sub(a): return a[idx]
z0,eps,r,kap,eta,phi,psi,delta,nu,sig,b0,be,bw,bf=map(sub,[z0,eps,r,kap,eta,phi,psi,delta,nu,sig,b0,be,bw,bf])
D=z0.shape[0]
K=np.exp(b0[:,None]+be[:,None]*elev+bw[:,None]*wet+bf[:,None]*forr)   # (D,m)
# reconstruct latent field M (D,m,T)
z=np.zeros((D,m,T)); z[:,:,0]=z0
for t in range(1,T):
    Mp=np.exp(z[:,:,t-1])
    surv=phi[:,None]*np.exp(-Mp/K)*Mp
    g=Mp/(Mp+nu[:,None]); src=(g*Mp)/deg
    immig=delta[:,None]*(src@W.T)
    mu=surv+immig+psi[:,None]+1e-9
    z[:,:,t]=np.log(mu)+sig[:,None]*eps[:,:,t-1]
M=np.exp(z)
totN=M.sum(1)                     # (D,T)
# ---------- Figure: latent field ----------
fig,ax=plt.subplots(1,3,figsize=(15,4))
tt=np.arange(T)
q=np.percentile(totN,[2.5,50,97.5],axis=0)
ax[0].fill_between(tt,q[0],q[2],alpha=.25,color='C0'); ax[0].plot(tt,q[1],'C0',label='posterior total $\\sum_i N$')
obsTot=np.nansum(np.where(obsC,YC,np.nan),axis=0)
ax[0].plot(tt,obsTot,'k.-',label='observed overnight total')
ax[0].set_yscale('symlog'); ax[0].set_xlabel('biweekly step'); ax[0].set_title('(a) total abundance'); ax[0].legend(fontsize=8)
peakt=int(np.argmax(q[1])); fieldpeak=M[:,:,peakt].mean(0)
ax[1].bar(np.arange(m),fieldpeak,color='C2'); ax[1].set_xlabel('plot'); ax[1].set_title('(b) mean field at peak step %d'%peakt)
for i in range(m): ax[2].plot(tt,M[:,i,:].mean(0),alpha=.7)
sc=np.where(obsC,YC,np.nan)
for i in range(m): ax[2].plot(tt,sc[i],'k.',ms=3)
ax[2].set_yscale('symlog'); ax[2].set_xlabel('biweekly step'); ax[2].set_title('(c) per-plot trajectories')
plt.tight_layout(); plt.savefig("mos_field.png",dpi=130); plt.close()
# ---------- Figure: posterior predictive ----------
nz=[]; mx=[]; dr=[]
for dd in range(D):
    YCrep=rng.poisson(rng.gamma(r[dd],np.maximum(M[dd],1e-6)/r[dd]))
    pi=eta[dd]+(1-eta[dd])*(1-np.exp(-kap[dd]*M[dd]))
    YPrep=(rng.random((m,T))<pi).astype(float)
    nz.append((YCrep[obsC]==0).mean()*obsC.sum()); mx.append(np.log10(max(YCrep[obsC].max(),1)))
    dr.append(YPrep[obsP].mean())
fig,ax=plt.subplots(1,3,figsize=(15,3.6))
ax[0].hist(nz,30,color='C0',alpha=.8); ax[0].axvline((YC[obsC]==0).sum(),color='r',lw=2); ax[0].set_title('(a) # zero counts')
ax[1].hist(mx,30,color='C0',alpha=.8); ax[1].axvline(np.log10(np.nanmax(YC)),color='r',lw=2); ax[1].set_title('(b) max count ($\\log_{10}$)')
ax[2].hist(dr,30,color='C0',alpha=.8); ax[2].axvline(np.nansum(YP)/obsP.sum(),color='r',lw=2); ax[2].set_title('(c) detection rate')
plt.tight_layout(); plt.savefig("mos_ppc.png",dpi=130); plt.close()
print("regenerated mos_field.png and mos_ppc.png; peak step",peakt,"obs detection rate %.2f"%(np.nansum(YP)/obsP.sum()))
