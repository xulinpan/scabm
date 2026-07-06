"""
Count-Valued Spatio-Temporal Agent-Based Model (SC-ABM)
=======================================================

Reference implementation of the generative model specified in
``SC_ABM_model_math_details.tex`` (X. Pan).

The latent abundance field N_{i,t} on a spatial lattice evolves as the sum of
three independent demographic components conditional on N_{t-1}:

    N_{i,t} = S_{i,t} + I_{i,t} + L_{i,t}

  * Density-regulated survival     S_{i,t} ~ Binomial(N_{i,t-1}, p_s)
        with  p_s = phi_i * exp(-N_{i,t-1}/K_i)                       (Ricker-type)
  * Neighbourhood immigration      I_{i,t} ~ Poisson(Lambda_{i,t})
        with  Lambda_{i,t} = sum_j  p_{j->i} * g(N_{j,t-1}) * N_{j,t-1}
        and   g(N;eta) = N/(N+eta)                                    (saturation)
  * Long-distance arrivals ("rain") L_{i,t} ~ Poisson(psi_i)

Spatial hierarchy (Sec. 5):
  * alpha ~ GP:  alpha = X beta + eps,  eps ~ N(0, sigma_a^2 C(theta_a))
                 C_ij = exp(-||s_i - s_j|| / theta_a)
  * log K_i     = x_i' beta_K + u_i,   u_i ~ N(0, sigma_K^2)
  * logit phi_i = x_i' beta_phi
  * log psi_i   = w_i' gamma_psi          (floored at eps0 > 0)

Dual-channel observation model (Sec. 6):
  * Count   channel: Y^C_{i,t} ~ NegBin(mean = rho * N, dispersion = r)
  * Binary  channel: Y^P_{i,t} ~ Bernoulli(pi),
                     pi = eta + (1-eta)(1 - exp(-kappa N))

Dispersal weights use the softmax-with-retention allocation (Eq. 8), which
strictly satisfies the unit dispersal-budget constraint  sum_i p_{j->i} <= 1.

The script also contains Monte-Carlo *verifications* of the paper's analytical
results:
  * one-step conditional moments (Prop. 3),
  * latent underdispersion  Var - E = -N p_s^2 <= 0  (Eq. 18),
  * the isolated-cell geometric-ergodicity stationary-mean bound
        E_pi[N] <= phi K / e + psi                                    (Thm. 12).

Author: generated for X. Pan's SC-ABM specification.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from scipy.spatial.distance import cdist


# =============================================================================
# 1. Model configuration
# =============================================================================
@dataclass
class ABMConfig:
    # --- spatial grid ---
    grid_rows: int = 15
    grid_cols: int = 15
    neighbourhood: str = "queen"      # "queen" (8-nn) or "rook" (4-nn)

    # --- temporal ---
    T: int = 60                       # number of time steps

    # --- covariate dimensions ---
    p: int = 2                        # habitat covariates x_i  (excl. intercept)
    q: int = 1                        # connectivity covariates w_i (excl. intercept)

    # --- hierarchy: habitat suitability GP  (Eq. 15) ---
    beta: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.4, -0.3]))
    sigma_alpha2: float = 0.25
    theta_alpha: float = 3.0

    # --- carrying capacity  log K = x'beta_K + u  (Eq. 16) ---
    beta_K: np.ndarray = field(default_factory=lambda: np.array([3.0, 0.5, 0.2]))
    sigma_K2: float = 0.10

    # --- survival  logit(phi) = x'beta_phi  (Eq. 17) ---
    beta_phi: np.ndarray = field(default_factory=lambda: np.array([0.8, 0.3, -0.2]))

    # --- long-distance rate  log psi = w'gamma_psi  (Eq. 10) ---
    gamma_psi: np.ndarray = field(default_factory=lambda: np.array([0.3, 0.25]))
    eps0: float = 0.05                # floor on psi_i (prevents absorbing 0 state)

    # --- dispersal ---
    eta: float = 0.5                  # density-scaling saturation constant
    disperse_scale: float = 1.0       # scale of softmax scores o_{ji}

    # --- observation: count channel (NegBin) ---
    rho: float = 0.6                  # detection fraction
    r: float = 5.0                    # NegBin dispersion (r -> inf => Poisson)

    # --- observation: binary channel (Bernoulli) ---
    kappa: float = 0.15               # per-individual detection intensity
    eta_fp: float = 0.02              # false-positive floor

    seed: int = 20260703


# =============================================================================
# 2. Grid + neighbourhood construction
# =============================================================================
def build_grid(cfg: ABMConfig):
    """Return cell coordinates and neighbour index lists for the lattice."""
    coords = np.array([(rr, cc)
                       for rr in range(cfg.grid_rows)
                       for cc in range(cfg.grid_cols)], dtype=float)
    m = coords.shape[0]

    if cfg.neighbourhood == "queen":
        offsets = [(-1, -1), (-1, 0), (-1, 1),
                   (0, -1),           (0, 1),
                   (1, -1),  (1, 0),  (1, 1)]
    elif cfg.neighbourhood == "rook":
        offsets = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    else:
        raise ValueError("neighbourhood must be 'queen' or 'rook'")

    def rc_to_idx(rr, cc):
        return rr * cfg.grid_cols + cc

    neighbours = [[] for _ in range(m)]
    for rr in range(cfg.grid_rows):
        for cc in range(cfg.grid_cols):
            i = rc_to_idx(rr, cc)
            for dr, dcc in offsets:
                nr, ncc = rr + dr, cc + dcc
                if 0 <= nr < cfg.grid_rows and 0 <= ncc < cfg.grid_cols:
                    neighbours[i].append(rc_to_idx(nr, ncc))
    return coords, neighbours


# =============================================================================
# 3. Hierarchy: draw covariates, GP suitability, and demographic parameters
# =============================================================================
def sample_hierarchy(cfg: ABMConfig, coords, rng):
    m = coords.shape[0]

    # covariates (with intercept column)
    X = np.column_stack([np.ones(m), rng.normal(size=(m, cfg.p))])          # m x (p+1)
    W = np.column_stack([np.ones(m), rng.normal(size=(m, cfg.q))])          # m x (q+1)

    # --- Gaussian-process latent suitability alpha (Eq. 15) ---
    D = cdist(coords, coords)                                               # Euclidean
    C = np.exp(-D / cfg.theta_alpha)
    Sigma = cfg.sigma_alpha2 * C + 1e-8 * np.eye(m)                         # jitter
    L_chol = np.linalg.cholesky(Sigma)
    alpha = X @ cfg.beta + L_chol @ rng.normal(size=m)

    # --- carrying capacity K_i (Eq. 16) ---
    u = rng.normal(scale=np.sqrt(cfg.sigma_K2), size=m)
    logK = X @ cfg.beta_K + u
    K = np.exp(logK)

    # --- survival baseline phi_i (Eq. 17); nudged by suitability alpha ---
    logit_phi = X @ cfg.beta_phi + 0.5 * (alpha - alpha.mean())
    phi = 1.0 / (1.0 + np.exp(-logit_phi))
    phi = np.clip(phi, 1e-4, 1 - 1e-4)

    # --- long-distance rate psi_i (Eq. 10) with floor eps0 ---
    logpsi = W @ cfg.gamma_psi
    psi = np.maximum(np.exp(logpsi), cfg.eps0)

    return dict(X=X, W=W, alpha=alpha, K=K, phi=phi, psi=psi)


# =============================================================================
# 4. Dispersal weights p_{j->i} via softmax-with-retention (Eq. 8)
# =============================================================================
def build_dispersal(cfg: ABMConfig, coords, neighbours, params, rng):
    """
    Returns p_out[j] = list of (target i, weight p_{j->i}) with
    sum_i p_{j->i} < 1 strictly (unit dispersal budget, Eq. 7).

    Scores o_{ji} favour targets that are close and more suitable:
        o_{ji} = disperse_scale * ( alpha_i - dist(j,i) )
    with retention score o_{j0} = 0.
    """
    m = coords.shape[0]
    alpha = params["alpha"]
    p_out = [[] for _ in range(m)]

    for j in range(m):
        nbrs = neighbours[j]
        if not nbrs:
            continue
        dists = np.linalg.norm(coords[nbrs] - coords[j], axis=1)
        scores = cfg.disperse_scale * (alpha[nbrs] - dists)
        # add mild noise so the network is non-degenerate
        scores = scores + 0.1 * rng.normal(size=len(nbrs))
        # softmax-with-retention: denominator includes exp(o_{j0}) = exp(0) = 1
        ex = np.exp(scores - scores.max())          # stabilise
        denom = np.exp(-scores.max()) + ex.sum()    # retention term exp(0)=1 -> exp(-max)
        weights = ex / denom
        for i, w in zip(nbrs, weights):
            p_out[j].append((i, float(w)))
    return p_out


def _incoming_map(m, p_out):
    """Invert p_out to p_in[i] = list of (source j, weight p_{j->i})."""
    p_in = [[] for _ in range(m)]
    for j, targets in enumerate(p_out):
        for (i, w) in targets:
            p_in[i].append((j, w))
    return p_in


# =============================================================================
# 5. Latent generative process  (Sec. 2)
# =============================================================================
def g_scale(N, eta):
    """Density-scaling g(N;eta) = N/(N+eta), with g(0)=0 (Eq. 6)."""
    return N / (N + eta)


def simulate_latent(cfg: ABMConfig, params, p_in, rng, N0=None):
    """
    Forward-simulate the latent abundance field N (shape T x m).
    """
    m = params["K"].shape[0]
    phi, K, psi = params["phi"], params["K"], params["psi"]

    N = np.zeros((cfg.T, m), dtype=np.int64)

    # --- initial field N_1: Poisson seeded near a fraction of carrying capacity ---
    if N0 is None:
        N[0] = rng.poisson(0.3 * K)
    else:
        N[0] = N0

    for t in range(1, cfg.T):
        Nprev = N[t - 1]

        # (a) density-regulated survival: Binomial(N_prev, p_s)
        p_s = phi * np.exp(-Nprev / K)
        p_s = np.clip(p_s, 0.0, 1.0)
        S = rng.binomial(Nprev, p_s)

        # (b) neighbourhood immigration rate Lambda_i  (Eq. 5)
        gN = g_scale(Nprev, cfg.eta) * Nprev                 # per-source emission
        Lambda = np.zeros(m)
        for i in range(m):
            acc = 0.0
            for (j, w) in p_in[i]:
                acc += w * gN[j]
            Lambda[i] = acc
        I = rng.poisson(Lambda)

        # (c) long-distance arrivals: Poisson(psi)
        Larr = rng.poisson(psi)

        N[t] = S + I + Larr

    return N


# =============================================================================
# 6. Dual-channel observation model  (Sec. 6)
# =============================================================================
def _nb_sample(mean, r, rng):
    """
    Sample NegBin in mean-dispersion parameterisation:
        Var = mean (1 + mean/r);  r -> inf => Poisson.
    Uses the Gamma-Poisson mixture: N ~ Poisson(Gamma(shape=r, scale=mean/r)).
    """
    mean = np.asarray(mean, dtype=float)
    out = np.zeros_like(mean)
    pos = mean > 0
    if np.any(pos):
        lam = rng.gamma(shape=r, scale=mean[pos] / r)
        out[pos] = rng.poisson(lam)
    return out.astype(np.int64)


def observe(cfg: ABMConfig, N, rng, count_frac=0.7, binary_frac=0.7):
    """
    Produce dual-channel observations. Each cell-time is observed in the count
    and/or binary channel with the given probabilities (missing => NaN / -1).
    """
    T, m = N.shape

    # count channel: Y^C ~ NegBin(rho N, r)
    Yc = _nb_sample(cfg.rho * N, cfg.r, rng).astype(float)
    obs_c = rng.random((T, m)) < count_frac
    Yc[~obs_c] = np.nan

    # binary channel: Y^P ~ Bernoulli(pi), pi = eta + (1-eta)(1 - e^{-kappa N})
    pi = cfg.eta_fp + (1 - cfg.eta_fp) * (1 - np.exp(-cfg.kappa * N))
    Yp = (rng.random((T, m)) < pi).astype(float)
    obs_p = rng.random((T, m)) < binary_frac
    Yp[~obs_p] = np.nan

    return dict(Yc=Yc, Yp=Yp, pi=pi, obs_c=obs_c, obs_p=obs_p)


# =============================================================================
# 7. Analytical-result verifications (Monte Carlo)
# =============================================================================
def verify_one_step_moments(cfg, n_prev=40, phi=0.7, K=50.0, Lambda=3.0,
                            psi=1.0, n_draws=400_000, rng=None):
    """
    Check Prop. 3: for the single-cell transition
        E[N_t | n]  = n p_s + Lambda + psi
        Var[N_t | n]= n p_s(1-p_s) + Lambda + psi
        Var - E     = -n p_s^2   (<= 0, underdispersion)
    """
    if rng is None:
        rng = np.random.default_rng(0)
    p_s = phi * np.exp(-n_prev / K)

    S = rng.binomial(n_prev, p_s, size=n_draws)
    M = rng.poisson(Lambda + psi, size=n_draws)     # I + L superposition
    Nt = S + M

    emp_mean, emp_var = Nt.mean(), Nt.var()
    thy_mean = n_prev * p_s + Lambda + psi
    thy_var = n_prev * p_s * (1 - p_s) + Lambda + psi
    return dict(p_s=p_s,
                emp_mean=emp_mean, thy_mean=thy_mean,
                emp_var=emp_var, thy_var=thy_var,
                emp_var_minus_mean=emp_var - emp_mean,
                thy_var_minus_mean=-n_prev * p_s ** 2)


def verify_isolated_stationary_bound(phi=0.7, K=50.0, psi=1.0,
                                     T=20_000, burn=2_000, rng=None):
    """
    Check Thm. 12: the isolated cell N_t = S_t + L_t is geometrically ergodic
    with stationary mean  E_pi[N] <= phi K / e + psi.
    """
    if rng is None:
        rng = np.random.default_rng(1)
    N = 0
    samples = np.empty(T)
    for t in range(T):
        p_s = phi * np.exp(-N / K)
        S = rng.binomial(N, np.clip(p_s, 0, 1))
        L = rng.poisson(psi)
        N = S + L
        samples[t] = N
    stat_mean = samples[burn:].mean()
    bound = phi * K / np.e + psi
    return dict(stat_mean=stat_mean, bound=bound, satisfied=stat_mean <= bound)


# =============================================================================
# 8. Driver
# =============================================================================
def run(cfg: ABMConfig | None = None, make_plots=True):
    if cfg is None:
        cfg = ABMConfig()
    rng = np.random.default_rng(cfg.seed)

    coords, neighbours = build_grid(cfg)
    m = coords.shape[0]
    params = sample_hierarchy(cfg, coords, rng)
    p_out = build_dispersal(cfg, coords, neighbours, params, rng)
    p_in = _incoming_map(m, p_out)

    # sanity: budget constraint sum_i p_{j->i} <= 1 for every source
    budgets = np.array([sum(w for (_, w) in p_out[j]) for j in range(m)])
    assert budgets.max() < 1.0 + 1e-9, "budget constraint violated"

    N = simulate_latent(cfg, params, p_in, rng)
    obs = observe(cfg, N, rng)

    # ---- console report ----
    print("=" * 70)
    print("SC-ABM simulation")
    print("=" * 70)
    print(f"grid            : {cfg.grid_rows} x {cfg.grid_cols}  = {m} cells "
          f"({cfg.neighbourhood} neighbourhood)")
    print(f"time steps      : {cfg.T}")
    print(f"dispersal budget: max_j sum_i p(j->i) = {budgets.max():.4f}  (< 1 OK)")
    print(f"K   range       : [{params['K'].min():.1f}, {params['K'].max():.1f}]")
    print(f"phi range       : [{params['phi'].min():.3f}, {params['phi'].max():.3f}]")
    print(f"psi range       : [{params['psi'].min():.3f}, {params['psi'].max():.3f}]")
    print(f"total abundance : t0 = {N[0].sum():,d}  ->  tT = {N[-1].sum():,d}")
    print(f"mean occupancy  : {(N[-1] > 0).mean():.3f}  (fraction non-empty)")

    print("\n--- Verify one-step conditional moments (Prop. 3) ---")
    mom = verify_one_step_moments(cfg, rng=np.random.default_rng(7))
    print(f"p_s              = {mom['p_s']:.4f}")
    print(f"E[N|n]   empirical={mom['emp_mean']:.4f}  theory={mom['thy_mean']:.4f}")
    print(f"Var[N|n] empirical={mom['emp_var']:.4f}  theory={mom['thy_var']:.4f}")
    print(f"Var - E  empirical={mom['emp_var_minus_mean']:.4f}  "
          f"theory={mom['thy_var_minus_mean']:.4f}  (<=0 underdispersed)")

    print("\n--- Verify isolated-cell stationary-mean bound (Thm. 12) ---")
    st = verify_isolated_stationary_bound(rng=np.random.default_rng(11))
    print(f"E_pi[N] empirical = {st['stat_mean']:.4f}")
    print(f"bound phiK/e + psi= {st['bound']:.4f}   satisfied = {st['satisfied']}")

    if make_plots:
        _plots(cfg, coords, params, N, obs)

    return dict(cfg=cfg, coords=coords, params=params, N=N, obs=obs,
                p_out=p_out, moments=mom, stationary=st)


# =============================================================================
# 9. Plots
# =============================================================================
def _plots(cfg, coords, params, N, obs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    R, Cc = cfg.grid_rows, cfg.grid_cols
    fig = plt.figure(figsize=(15, 10))

    # (1) total abundance over time
    ax = fig.add_subplot(2, 3, 1)
    ax.plot(N.sum(axis=1), lw=2, color="C0")
    ax.set_title("Total latent abundance  $\\sum_i N_{i,t}$")
    ax.set_xlabel("time $t$"); ax.set_ylabel("total abundance")
    ax.grid(alpha=0.3)

    # (2) per-cell trajectories (a random subset)
    ax = fig.add_subplot(2, 3, 2)
    rng = np.random.default_rng(0)
    for i in rng.choice(N.shape[1], size=12, replace=False):
        ax.plot(N[:, i], lw=0.9, alpha=0.7)
    ax.set_title("Sample per-cell trajectories $N_{i,t}$")
    ax.set_xlabel("time $t$"); ax.set_ylabel("abundance")
    ax.grid(alpha=0.3)

    # (3) carrying capacity field K_i
    ax = fig.add_subplot(2, 3, 3)
    im = ax.imshow(params["K"].reshape(R, Cc), origin="lower", cmap="viridis")
    ax.set_title("Carrying capacity $K_i$"); fig.colorbar(im, ax=ax, shrink=0.8)

    # (4) latent field snapshot at final time
    ax = fig.add_subplot(2, 3, 4)
    im = ax.imshow(N[-1].reshape(R, Cc), origin="lower", cmap="magma")
    ax.set_title(f"Latent field $N_{{i,T}}$  (t={cfg.T-1})")
    fig.colorbar(im, ax=ax, shrink=0.8)

    # (5) count-channel observation snapshot
    ax = fig.add_subplot(2, 3, 5)
    Yc_last = obs["Yc"][-1].reshape(R, Cc)
    im = ax.imshow(np.ma.masked_invalid(Yc_last), origin="lower", cmap="cividis")
    ax.set_title("Count obs $Y^{(C)}_{i,T}$  (NB, missing=grey)")
    ax.set_facecolor("lightgrey")
    fig.colorbar(im, ax=ax, shrink=0.8)

    # (6) binary detection probability field pi
    ax = fig.add_subplot(2, 3, 6)
    im = ax.imshow(obs["pi"][-1].reshape(R, Cc), origin="lower",
                   cmap="plasma", vmin=0, vmax=1)
    ax.set_title("Detection prob $\\pi_{i,T}=\\eta+(1-\\eta)(1-e^{-\\kappa N})$")
    fig.colorbar(im, ax=ax, shrink=0.8)

    fig.suptitle("Count-Valued Spatio-Temporal ABM — simulation", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = "sc_abm_simulation.png"
    fig.savefig(out, dpi=130)
    print(f"\n[figure saved to {out}]")

    # second figure: moment-matching diagnostics
    fig2, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # underdispersion sweep: Var - E vs n_prev
    ns = np.arange(0, 120, 4)
    vmm_emp, vmm_thy = [], []
    r2 = np.random.default_rng(3)
    for n0 in ns:
        mom = verify_one_step_moments(cfg, n_prev=n0, n_draws=120_000, rng=r2)
        vmm_emp.append(mom["emp_var_minus_mean"])
        vmm_thy.append(mom["thy_var_minus_mean"])
    axes[0].plot(ns, vmm_thy, "k-", lw=2, label="theory $-N p_s^2$")
    axes[0].plot(ns, vmm_emp, "o", ms=4, color="C3", label="Monte Carlo")
    axes[0].axhline(0, color="grey", lw=0.8)
    axes[0].set_title("Latent underdispersion (Prop. 3):  Var $-$ E $\\leq 0$")
    axes[0].set_xlabel("$N_{i,t-1}$"); axes[0].set_ylabel("Var $-$ E")
    axes[0].legend(); axes[0].grid(alpha=0.3)

    # stationary distribution of isolated cell vs. mean bound
    r3 = np.random.default_rng(5)
    N1 = 0; samp = []
    for _ in range(40_000):
        ps = 0.7 * np.exp(-N1 / 50.0)
        N1 = r3.binomial(N1, np.clip(ps, 0, 1)) + r3.poisson(1.0)
        samp.append(N1)
    samp = np.array(samp[2000:])
    axes[1].hist(samp, bins=40, density=True, color="C0", alpha=0.7)
    bound = 0.7 * 50.0 / np.e + 1.0
    axes[1].axvline(samp.mean(), color="C3", lw=2,
                    label=f"$E_\\pi[N]$={samp.mean():.1f}")
    axes[1].axvline(bound, color="k", lw=2, ls="--",
                    label=f"bound $\\phi K/e+\\psi$={bound:.1f}")
    axes[1].set_title("Isolated-cell stationary law (Thm. 12)")
    axes[1].set_xlabel("$N$"); axes[1].set_ylabel("density")
    axes[1].legend(); axes[1].grid(alpha=0.3)

    fig2.tight_layout()
    out2 = "sc_abm_diagnostics.png"
    fig2.savefig(out2, dpi=130)
    print(f"[figure saved to {out2}]")


if __name__ == "__main__":
    run()
