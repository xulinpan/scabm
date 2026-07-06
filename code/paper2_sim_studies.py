"""
Paper 2 (JABES) revision simulation studies — self-contained, exact-HMM-filter
Laplace analyses addressing editor points 2.3-2.5. Checkpointed to results file.

  A (2.3 joint process): observation-vector coverage with (phi,psi,K) FIXED at
     truth vs JOINTLY estimated.
  B (2.4 low abundance):  coverage in a low/high-zero regime (mean~1.5) vs moderate.
  C (2.5 eta prior):      eta identification (true eta=0.15) under informative,
     diffuse, and mis-centred priors.

Estimator: exact truncated forward (HMM) filter marginalising the integer latent
field per cell (isolated-cell = the conditional regime of Prop 3), MAP + Laplace
(observed-information) 95% Wald intervals. Deterministic given seeds.

Usage:  python paper2_sim_studies.py            # runs a time-boxed chunk, checkpoints
        python paper2_sim_studies.py summarize  # prints coverage table
"""
import sys, os, json, time, numpy as np
from math import lgamma
from scipy.optimize import minimize
from scipy.stats import binom, poisson
from scipy.special import gammaln as _lg

NMAX = 120
STATES = np.arange(NMAX + 1)
RESULTS = "paper2_sim_results.jsonl"
R, M, T, BUDGET = 20, 35, 18, 38.0

def expit(x): return 1.0/(1.0+np.exp(-x))

def transition_matrix(phi, K, psi):
    Tm = np.zeros((NMAX+1, NMAX+1)); pois = poisson.pmf(STATES, psi)
    for n in range(NMAX+1):
        ps = phi*np.exp(-n/K); b = binom.pmf(np.arange(n+1), n, ps)
        conv = np.convolve(b, pois)[:NMAX+1]
        row = np.zeros(NMAX+1); row[:len(conv)] = conv
        row[NMAX] += max(0.0, 1.0-row.sum()); Tm[n] = row
    return Tm
def stationary_p0(Tm, it=60):
    p = np.zeros(NMAX+1); p[1] = 1.0
    for _ in range(it): p = p @ Tm
    return p
_CACHE = {}
def trans_p0(phi, K, psi):
    key = (round(phi,5), round(K,4), round(psi,5)); h = _CACHE.get(key)
    if h is None:
        Tm = transition_matrix(phi,K,psi); h = (Tm, stationary_p0(Tm))
        if len(_CACHE) < 30000: _CACHE[key] = h
    return h

def simulate(th, m, T, seed):
    rr = np.random.default_rng(seed)
    rho, r, kap, eta = expit(th['mu_rho']), th['r'], th['kappa'], th['eta']
    phi, psi, K = th['phi'], th['psi'], th['K']
    N = np.zeros((m, T+1), int); N[:,0] = rr.poisson(2, m)
    for t in range(1, T+1):
        ps = phi*np.exp(-N[:,t-1]/K)
        N[:,t] = rr.binomial(N[:,t-1], ps) + rr.poisson(psi, m)
    N = N[:,1:]; mean = rho*N + 1e-9
    YC = np.minimum(rr.poisson(rr.gamma(r, mean/r)), NMAX*40)
    pi = eta + (1-eta)*(1-np.exp(-kap*N)); YP = (rr.random((m,T)) < pi).astype(int)
    return YC, YP, N

def count_obs(YC, mu_rho, r):
    rho = expit(mu_rho); mean = rho*STATES + 1e-9
    logC = _lg(YC[...,None]+r) - lgamma(r) - _lg(YC[...,None]+1.0)
    return np.exp(logC + r*np.log(r/(r+mean))[None,None,:] + YC[...,None]*np.log(mean/(r+mean))[None,None,:])
def binary_obs(YP, kappa, eta):
    pi = np.clip(eta + (1-eta)*(1-np.exp(-kappa*STATES)), 1e-9, 1-1e-9)
    return np.where(YP[...,None]==1, pi[None,None,:], 1-pi[None,None,:])
def forward(Tm, Obs, p0):
    A = p0[None,:]*Obs[:,0,:]; c = A.sum(1); ll = np.log(c).sum(); A /= c[:,None]
    for t in range(1, Obs.shape[1]):
        A = (A@Tm)*Obs[:,t,:]; c = A.sum(1); ll += np.log(c).sum(); A /= c[:,None]
    return ll

def log_prior(th, ep):
    lp = -0.5*(th['mu_rho']/1.5)**2 + (2-1)*np.log(th['r']) - 0.5*th['r'] \
         - 0.5*(th['kappa']/0.1)**2
    a,b = ep; lp += (a-1)*np.log(th['eta']) + (b-1)*np.log(1-th['eta'])
    lp += -0.5*((np.log(th['phi'])-np.log(0.5))/1.0)**2
    lp += -0.5*((np.log(th['psi'])-np.log(4.0))/1.0)**2
    lp += -0.5*((np.log(th['K'])-np.log(40.0))/1.0)**2
    return lp

FREE_ALL = ['mu_rho','r','kappa','eta','phi','psi','K']
BND = {'mu_rho':(-6,6),'r':(np.log(0.1),np.log(50)),'kappa':(np.log(1e-3),np.log(5)),
       'eta':(-9,1),'phi':(-6,6),'psi':(np.log(1e-2),np.log(50)),'K':(np.log(2),np.log(400))}
def pack(th, free):
    v=[]
    for k in free:
        x=th[k]
        v.append(x if k=='mu_rho' else np.log(x/(1-x)) if k in ('eta','phi') else np.log(x))
    return np.array(v)
def unpack(v, free, base):
    th=dict(base)
    for k,x in zip(free,v):
        th[k]= x if k=='mu_rho' else expit(x) if k in ('eta','phi') else np.exp(x)
    return th
def nlp(v, free, base, YC, YP, ep):
    th=unpack(v,free,base)
    if th['r']<=0 or th['kappa']<=0 or th['psi']<=0 or th['K']<=0: return 1e12
    Tm,p0=trans_p0(th['phi'],th['K'],th['psi'])
    Obs=count_obs(YC,th['mu_rho'],th['r'])*binary_obs(YP,th['kappa'],th['eta'])
    return -(forward(Tm,Obs,p0)+log_prior(th,ep))
def laplace(free, base, truth, YC, YP, ep):
    v0=pack(truth,free); bnds=[BND[k] for k in free]
    res=minimize(nlp,v0,args=(free,base,YC,YP,ep),method='L-BFGS-B',bounds=bnds,
                 options={'maxiter':200,'ftol':1e-7})
    v=res.x; n=len(v); H=np.zeros((n,n)); e=1e-3
    for i in range(n):
        for j in range(i,n):
            ei=np.zeros(n);ei[i]=e; ej=np.zeros(n);ej[j]=e
            H[i,j]=H[j,i]=(nlp(v+ei+ej,free,base,YC,YP,ep)-nlp(v+ei-ej,free,base,YC,YP,ep)
                          -nlp(v-ei+ej,free,base,YC,YP,ep)+nlp(v-ei-ej,free,base,YC,YP,ep))/(4*e*e)
    try: sd=np.sqrt(np.clip(np.diag(np.linalg.inv(H)),0,None))
    except np.linalg.LinAlgError: sd=np.full(n,np.nan)
    return v, sd
def natural(k, xh, xs):
    lo,hi = xh-1.96*xs, xh+1.96*xs
    f = (lambda z:z) if k=='mu_rho' else (expit if k in ('eta','phi') else np.exp)
    return f(xh), f(lo), f(hi)

TM = dict(mu_rho=-1.0, r=3.0, kappa=0.08, eta=0.02, phi=0.6, psi=4.0, K=40.0)
TL = dict(mu_rho=-1.0, r=3.0, kappa=0.20, eta=0.02, phi=0.5, psi=1.0, K=6.0)
TF = dict(mu_rho=-1.0, r=3.0, kappa=0.08, eta=0.15, phi=0.6, psi=4.0, K=40.0)
O4 = ['mu_rho','r','kappa','eta']
JOBS = [
    ("A1_fixed", TM, O4,       (1,49), ['mu_rho','kappa'], 100),
    ("A2_joint", TM, FREE_ALL, (1,49), ['mu_rho','kappa'], 100),
    ("B1_low",   TL, O4,       (1,49), ['mu_rho','kappa'], 200),
    ("B2_mod",   TM, O4,       (1,49), ['mu_rho','kappa'], 200),
    ("C_inform", TF, O4, (3,17), ['eta','kappa'], 300),
    ("C_diffuse",TF, O4, (1,1),  ['eta','kappa'], 300),
    ("C_miscen", TF, O4, (1,49), ['eta','kappa'], 300),
]
def done_set():
    s=set()
    if os.path.exists(RESULTS):
        for l in open(RESULTS):
            try: d=json.loads(l); s.add((d['job'],d['rep']))
            except: pass
    return s
def run():
    done=done_set(); t0=time.time(); n=0
    with open(RESULTS,"a") as f:
        for (name,truth,free,ep,targets,seed0) in JOBS:
            for rep in range(R):
                if (name,rep) in done: continue
                if time.time()-t0>BUDGET: print(f"budget; +{n} fits"); return
                YC,YP,N=simulate(truth,M,T,seed0+rep)
                try: v,sd=laplace(free,dict(truth),truth,YC,YP,ep)
                except Exception as e: print("FAIL",name,rep,repr(e)); continue
                rec={'job':name,'rep':rep,'meanN':float(N.mean()),'zf':float((YC==0).mean())}
                for k in targets:
                    i=free.index(k); est,lo,hi=natural(k,v[i],sd[i])
                    rec[k]={'cover':bool(np.isfinite(lo) and lo<=truth[k]<=hi),
                            'w':float(hi-lo) if np.isfinite(hi) else None}
                f.write(json.dumps(rec)+"\n"); f.flush(); n+=1
    print(f"ALL DONE (+{n} fits)")
def summarize():
    from collections import defaultdict
    recs=[json.loads(l) for l in open(RESULTS)]; by=defaultdict(list)
    for r in recs: by[r['job']].append(r)
    print(f"{'job':10s}{'n':>4}{'meanN':>7}{'zero%':>7}   targets")
    for (name,truth,free,ep,targets,seed0) in JOBS:
        rs=by.get(name,[])
        if not rs: print(f"{name:10s} (none)"); continue
        n=len(rs); mN=np.mean([x['meanN'] for x in rs]); zf=100*np.mean([x['zf'] for x in rs])
        parts=[]
        for k in targets:
            cs=[x[k]['cover'] for x in rs if k in x]; ws=[x[k]['w'] for x in rs if k in x and x[k]['w']]
            c=np.mean(cs); mcse=(c*(1-c)/max(len(cs),1))**0.5
            parts.append(f"{k} cov {c:.2f}(±{mcse:.2f}) w {np.mean(ws):.3f}" if ws else f"{k} cov {c:.2f}")
        print(f"{name:10s}{n:>4}{mN:>7.1f}{zf:>7.1f}   " + " | ".join(parts))
if __name__=="__main__":
    summarize() if (len(sys.argv)>1 and sys.argv[1]=="summarize") else run()
