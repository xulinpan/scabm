"""Run NEON NUTS chains one at a time (per invocation) to fit the per-call time
budget; each call runs the next unfinished chain and saves it. Then combine.py
computes split-R-hat and summaries."""
import sys, os, glob, time, numpy as np, jax, numpyro
numpyro.set_host_device_count(1)
import mos_numpyro_fit as MF
from numpyro.infer import MCMC, NUTS
NW, NS, NCHAIN = 400, 400, 4
done = sorted(glob.glob("mos_chain_*.npz"))
nxt = len(done)
if nxt >= NCHAIN:
    print("ALL CHAINS DONE", nxt); sys.exit()
t0 = time.time()
mcmc = MCMC(NUTS(MF.model, target_accept_prob=0.9, max_tree_depth=10),
            num_warmup=NW, num_samples=NS, num_chains=1, progress_bar=False)
mcmc.run(jax.random.PRNGKey(100 + nxt))
post = mcmc.get_samples()
np.savez(f"mos_chain_{nxt}.npz", **{k: np.asarray(v) for k, v in post.items()})
print(f"chain {nxt} done in {time.time()-t0:.1f}s ({nxt+1}/{NCHAIN})")
