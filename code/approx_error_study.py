"""
Approximation-error characterization for SC-ABM-NKD.

Two approximations are used at inference time:
  (L1) log-normal relaxation of integer latent counts  (main text Eq. 21)
  (L2) moment-matched Gaussian transition density        (main text Eq. 23)

Both replace the exact one-step transition
      X = S + A,   S ~ Binomial(n, p_s),   A ~ Poisson(a),   a = Lambda + psi,
      mu  = n p_s + a,   sigma^2 = n p_s (1 - p_s) + a,
by a continuous, moment-matched surrogate. This script quantifies the
resulting discrepancy as a function of the conditional mean count mu, and
shows that the per-cell error is governed by local mean count, not grid size.

Everything here is a self-contained numerical computation from the model
definition; it does not use the SGVB estimator, so it is fully reproducible.
"""
import numpy as np
from scipy.stats import binom, poisson, norm, lognorm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(20260701)


# ----------------------------------------------------------------------
# Exact transition pmf and moment-matched continuous surrogates
# ----------------------------------------------------------------------
def exact_pmf(n, ps, a, K):
    s = np.arange(0, n + 1)
    binom_pmf = binom.pmf(s, n, ps)
    kA = np.arange(0, K + 1)
    pois_pmf = poisson.pmf(kA, a)
    conv = np.convolve(binom_pmf, pois_pmf)
    p = np.zeros(K + 1)
    m = min(len(conv), K + 1)
    p[:m] = conv[:m]
    return p


def gaussian_discrete(mu, sig, K):
    k = np.arange(0, K + 1)
    p = norm.cdf((k + 0.5 - mu) / sig) - norm.cdf((k - 0.5 - mu) / sig)
    p[0] += norm.cdf((-0.5 - mu) / sig)  # fold sub-zero mass onto 0
    return p


def lognormal_discrete(mu, var, K):
    s2 = np.log(1.0 + var / mu ** 2)
    mln = np.log(mu) - 0.5 * s2
    s = np.sqrt(s2)
    cdf = lambda x: lognorm.cdf(x, s=s, scale=np.exp(mln))
    k = np.arange(0, K + 1)
    lo = np.maximum(k - 0.5, 0.0)
    p = cdf(k + 0.5) - cdf(lo)
    return p


def tv(p, q):
    n = max(len(p), len(q))
    P, Q = np.zeros(n), np.zeros(n)
    P[: len(p)] = p
    Q[: len(q)] = q
    return 0.5 * np.abs(P - Q).sum()


def kl(p, q, eps=1e-12):
    n = max(len(p), len(q))
    P, Q = np.zeros(n), np.zeros(n)
    P[: len(p)] = p
    Q[: len(q)] = q
    mask = P > eps
    return float(np.sum(P[mask] * np.log(P[mask] / np.clip(Q[mask], eps, None))))


# ----------------------------------------------------------------------
# Part A/C: error vs conditional mean count mu, disentangling L1 and L2
# ----------------------------------------------------------------------
def error_vs_mu(ps_values=(0.3, 0.6, 0.9), a=2.0, n_grid=None):
    if n_grid is None:
        n_grid = np.unique(np.round(np.linspace(1, 120, 60)).astype(int))
    rows = []
    for ps in ps_values:
        for n in n_grid:
            mu = n * ps + a
            var = n * ps * (1 - ps) + a
            sig = np.sqrt(var)
            K = int(mu + 12 * sig + 40)
            pe = exact_pmf(n, ps, a, K)
            pg = gaussian_discrete(mu, sig, K)
            pl = lognormal_discrete(mu, var, K)
            rows.append(dict(ps=ps, n=n, mu=mu, sigma=sig,
                             tv_gauss=tv(pe, pg), tv_lnorm=tv(pe, pl),
                             kl_gauss=kl(pe, pg), kl_lnorm=kl(pe, pl),
                             p0_exact=pe[0], p0_gauss=pg[0], p0_lnorm=pl[0]))
    return rows


# ----------------------------------------------------------------------
# Part B: per-cell error is governed by local mean, not grid size.
# Forward-simulate the exact process on grids of increasing size and
# record the realized per-cell TV of the Gaussian transition surrogate.
# ----------------------------------------------------------------------
def simple_kernel_weights(coords, decay=0.9):
    # fixed, isotropic Queen-neighbourhood weights (kernel form is irrelevant
    # to the transition-approximation error, which depends only on n, p_s, a)
    m = len(coords)
    W = np.zeros((m, m))
    for i in range(m):
        d = np.abs(coords - coords[i]).sum(1)
        nb = (d == 1)
        W[i, nb] = decay / max(nb.sum(), 1)
    return W


def grid_forward_study(sizes=(25, 100, 225, 400), T=8, phi=0.75, Kcap=14.0,
                       psi=1.5, eta=8.0):
    out = []
    for m in sizes:
        side = int(round(np.sqrt(m)))
        coords = np.array([(i, j) for i in range(side) for j in range(side)])
        m = len(coords)
        W = simple_kernel_weights(coords)
        N = rng.poisson(9.6, size=m).astype(int)  # app abundance regime
        tvs = []
        for t in range(T):
            ps = phi * np.exp(-N / Kcap)
            g = N / (N + eta)
            Lam = W @ (g * N)
            a = Lam + psi
            # record per-cell transition-approx TV at current state
            for i in range(m):
                n_i = int(N[i])
                mu = n_i * ps[i] + a[i]
                var = n_i * ps[i] * (1 - ps[i]) + a[i]
                sig = np.sqrt(max(var, 1e-9))
                K = int(mu + 12 * sig + 40)
                pe = exact_pmf(n_i, ps[i], a[i], K)
                pg = gaussian_discrete(mu, sig, K)
                tvs.append(tv(pe, pg))
            # advance exactly
            S = rng.binomial(N, ps)
            A = rng.poisson(a)
            N = S + A
        tvs = np.array(tvs)
        out.append(dict(m=m, mean_tv=tvs.mean(), se_tv=tvs.std() / np.sqrt(len(tvs)),
                        med_tv=np.median(tvs)))
    return out


def main():
    rowsA = error_vs_mu()
    gridB = grid_forward_study()

    # ---- console summary table (Part A/C at representative mu) ----
    print("mu     ps    TV_gauss   TV_lnorm   KL_gauss   KL_lnorm")
    for r in rowsA:
        if r["ps"] == 0.6 and int(round(r["mu"])) in (3, 6, 10, 20, 40, 60):
            print(f"{r['mu']:5.1f} {r['ps']:.1f}  {r['tv_gauss']:.4f}    "
                  f"{r['tv_lnorm']:.4f}    {r['kl_gauss']:.4f}    {r['kl_lnorm']:.4f}")
    print("\ngrid m   mean per-cell TV (Gaussian)   SE")
    for r in gridB:
        print(f"{r['m']:5d}      {r['mean_tv']:.4f}                 {r['se_tv']:.4f}")

    # ---- figure ----
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for ps, c in zip((0.3, 0.6, 0.9), ("#1b9e77", "#d95f02", "#7570b3")):
        sub = [r for r in rowsA if r["ps"] == ps]
        mu = [r["mu"] for r in sub]
        ax[0].plot(mu, [r["tv_gauss"] for r in sub], "-", color=c,
                   label=f"Gaussian, $p_s$={ps}")
        ax[0].plot(mu, [r["tv_lnorm"] for r in sub], "--", color=c,
                   label=f"log-normal, $p_s$={ps}")
    # reference mu^{-1/2} slope
    mm = np.linspace(3, 60, 50)
    ref = 0.20 * mm ** -0.5
    ax[0].plot(mm, ref, ":", color="gray", lw=1.4, label=r"$\propto \mu^{-1/2}$")
    ax[0].set_xlabel(r"conditional mean count $\mu$")
    ax[0].set_ylabel("total-variation distance to exact transition")
    ax[0].set_title("(a) Transition-approximation error vs. mean count")
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].legend(fontsize=7, ncol=2)
    ax[0].grid(alpha=0.3, which="both")

    ms = [r["m"] for r in gridB]
    mtv = [r["mean_tv"] for r in gridB]
    setv = [r["se_tv"] for r in gridB]
    ax[1].errorbar(ms, mtv, yerr=setv, fmt="o-", color="#d95f02", capsize=4)
    ax[1].set_xlabel("grid size $m$ (cells)")
    ax[1].set_ylabel("mean per-cell TV (Gaussian)")
    ax[1].set_title("(b) Per-cell error is flat in grid size\n(app. abundance regime)")
    ax[1].set_ylim(0, max(mtv) * 1.6)
    ax[1].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig("figs/figS_approx_error.pdf", bbox_inches="tight")
    print("\nsaved figs/figS_approx_error.pdf")

    # emit LaTeX-ready numbers
    def at(mu_target, ps=0.6, key="tv_gauss"):
        sub = [r for r in rowsA if r["ps"] == ps]
        r = min(sub, key=lambda x: abs(x["mu"] - mu_target))
        return r[key], r["mu"]
    print("\n-- numbers for text --")
    for mt in (3, 10, 30):
        vg, mu = at(mt, key="tv_gauss")
        vl, _ = at(mt, key="tv_lnorm")
        print(f"mu~{mu:.1f}: TV_gauss={vg:.3f}, TV_lnorm={vl:.3f}")


if __name__ == "__main__":
    main()
