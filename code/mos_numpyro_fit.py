"""
Definitive dual-channel spatial fit for the CORRECTED (species-specific) NEON
UNDE 2024 Coquillettidia perturbans data, in NumPyro (NUTS).

Run on a machine with numpyro + jax installed:
    pip install numpyro
    python mos_numpyro_fit.py

Reads mos_arrays.npz (produced by build_mos_arrays.py). Prints posterior
summaries with split-R-hat and ESS, and saves posterior draws to mos_numpyro_posterior.npz.

Model (relaxed abundance M_{i,t}=exp(z_{i,t}) > 0):
  transition   z_{i,t} = log mu_{i,t} + sigma_p * eps_{i,t},  eps ~ N(0,1)   (non-centred)
     mu_{i,t} = phi e^{-M_{i,t-1}/K_i} M_{i,t-1}
                + delta * sum_{j~i} g(M_{j,t-1};nu) M_{j,t-1}/deg_j + psi
     log K_i  = b0 + b_e elev_i + b_w wet_i + b_f for_i,   g(x;nu)=x/(x+nu)
  count   Y^C_{i,t} ~ NegBinomial2(mean = M_{i,t}, concentration = r)   (rho = 1, scale anchor)
  binary  Y^P_{i,t} ~ Bernoulli( eta + (1-eta)(1-e^{-kappa M_{i,t}}) )

Scale note: rho is fixed to 1 to anchor the abundance scale to the overnight
count channel (Section 6). For a data-driven scale, free rho and place an
informative prior on mu_rho (e.g. Normal(0,0.5)); the absolute scale is only
weakly identified otherwise (see revision notes).
"""
import numpy as np, jax, jax.numpy as jnp
import numpyro, numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS
import os

# Set FREE_RHO=1 to estimate rho under an informative prior instead of fixing rho=1
FREE_RHO = os.environ.get('FREE_RHO', '0') == '1'

d = np.load("mos_arrays.npz", allow_pickle=True)
YC = d['YC'].astype(float); YP = d['YP'].astype(float)
m, T = YC.shape
elev = np.nan_to_num((d['elev'] - np.nanmean(d['elev'])) / (np.nanstd(d['elev']) + 1e-9))
wet = d['wetland'].astype(float); forr = d['forest'].astype(float)
W = jnp.asarray(d['W'].astype(float)); deg = jnp.asarray(np.maximum(d['W'].sum(1), 1.0))
obsC = jnp.asarray(np.isfinite(YC)); obsP = jnp.asarray(np.isfinite(YP))
YCj = jnp.asarray(np.where(np.isfinite(YC), YC, 0.0))
YPj = jnp.asarray(np.where(np.isfinite(YP), YP, 0.0))
Y0 = jnp.asarray(np.where(np.isfinite(YC[:, 0]), np.nan_to_num(YC[:, 0]), 1.0))
elev, wet, forr = map(jnp.asarray, (elev, wet, forr))

def model():
    r     = numpyro.sample("r",     dist.Gamma(2.0, 0.5))
    kappa = numpyro.sample("kappa", dist.HalfNormal(0.1))
    eta   = numpyro.sample("eta",   dist.Beta(1.0, 49.0))
    phi   = numpyro.sample("phi",   dist.Beta(2.0, 2.0))
    psi   = numpyro.sample("psi",   dist.LogNormal(jnp.log(4.0), 1.0))
    delta = numpyro.sample("delta", dist.Beta(2.0, 2.0))
    nu    = numpyro.sample("nu",    dist.LogNormal(jnp.log(40.0), 1.0))
    sigp  = numpyro.sample("sigma_p", dist.HalfNormal(1.0))
    b0 = numpyro.sample("b0", dist.Normal(jnp.log(jnp.maximum(jnp.nanmean(YCj), 1.0)), 3.0))
    be = numpyro.sample("b_elev", dist.Normal(0.0, 2.0))
    bw = numpyro.sample("b_wet",  dist.Normal(0.0, 2.0))
    bf = numpyro.sample("b_for",  dist.Normal(0.0, 2.0))
    if FREE_RHO:
        mu_rho = numpyro.sample("mu_rho", dist.Normal(0.0, 0.5))  # informative scale prior
        rho = jax.nn.sigmoid(mu_rho)
    else:
        rho = 1.0
    K = jnp.exp(b0 + be*elev + bw*wet + bf*forr)

    z0  = numpyro.sample("z0",  dist.Normal(jnp.log(Y0 + 2.0), 2.0))
    eps = numpyro.sample("eps", dist.Normal(0, 1).expand([m, T-1]).to_event(2))

    def stepfn(z_prev, e):
        M_prev = jnp.exp(z_prev)
        surv = phi * jnp.exp(-M_prev / K) * M_prev
        g = M_prev / (M_prev + nu)
        immig = delta * (W @ ((g * M_prev) / deg))
        mu = surv + immig + psi + 1e-9
        z_t = jnp.log(mu) + sigp * e
        return z_t, z_t
    _, zs = jax.lax.scan(stepfn, z0, eps.T)      # zs: (T-1, m)
    z = jnp.concatenate([z0[None, :], zs], axis=0).T   # (m,T)
    M = jnp.exp(z)

    with numpyro.handlers.mask(mask=obsC):
        numpyro.sample("YC", dist.NegativeBinomial2(mean=rho * M + 1e-9, concentration=r), obs=YCj)
    pi = jnp.clip(eta + (1 - eta) * (1 - jnp.exp(-kappa * M)), 1e-6, 1 - 1e-6)
    with numpyro.handlers.mask(mask=obsP):
        numpyro.sample("YP", dist.Bernoulli(probs=pi), obs=YPj)
    numpyro.deterministic("total_N", M.sum(0))

if __name__ == "__main__":
    numpyro.set_host_device_count(4)
    kernel = NUTS(model, target_accept_prob=0.9, max_tree_depth=12)
    mcmc = MCMC(kernel, num_warmup=1500, num_samples=1500, num_chains=4, progress_bar=True)
    mcmc.run(jax.random.PRNGKey(0))
    mcmc.print_summary()
    post = mcmc.get_samples(group_by_chain=False)
    np.savez("mos_numpyro_posterior.npz", **{k: np.asarray(v) for k, v in post.items()})
    print("saved mos_numpyro_posterior.npz")
