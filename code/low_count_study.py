"""
Low-count robustness study for SC-ABM-NKD (addresses reviewer points 3.1/3.3:
extend the approximation validation to the low-abundance regime, mean count
~1-2, where both continuous relaxations are worst, and report agreement on the
LATENT FIELD, not only on scalar parameters).

Part 1 (exact, self-contained). One-step transition-approximation error in the
low-count regime. For an isolated cell we compare the EXACT Binomial-Poisson
transition against (i) the moment-matched Gaussian transition (Eq. 23) and
(ii) the moment-matched log-normal relaxation (Eq. 21) via total variation,
KL, the probability mass at zero P(N=0) (the ecologically important quantity
the log-normal cannot represent), and the calibration of each surrogate's
nominal central-90% interval measured against the exact transition.

Part 2 (particle filter). Latent-field recovery on a small COUPLED grid. A
bootstrap particle filter with the EXACT transition is the gold standard; a
second filter identical except that it propagates with the moment-matched
Gaussian transition isolates the effect of the transition relaxation on the
latent-field posterior. We report posterior-mean RMSE against the true field
and empirical 90% credible-interval coverage, in a low-abundance regime
(mean count ~1.5) and in the application regime (mean count ~9.6).
Self-contained (numpy + scipy).
"""
import numpy as np
from math import lgamma
from scipy.stats import binom, poisson, norm, lognorm

rng_global = np.random.default_rng(20260702)

# ----------------------------------------------------------------------
# Part 1: exact transition vs surrogates in the low-count regime
# ----------------------------------------------------------------------
def exact_pmf(n, ps, a, K):
    s = np.arange(0, n + 1)
    conv = np.convolve(binom.pmf(s, n, ps), poisson.pmf(np.arange(0, K + 1), a))
    p = np.zeros(K + 1); m = min(len(conv), K + 1); p[:m] = conv[:m]
    return p

def gaussian_discrete(mu, sig, K):
    k = np.arange(0, K + 1)
    p = norm.cdf((k + 0.5 - mu) / sig) - norm.cdf((k - 0.5 - mu) / sig)
    p[0] += norm.cdf((-0.5 - mu) / sig)          # fold sub-zero mass onto 0
    return p

def lognormal_discrete(mu, var, K):
    s2 = np.log(1.0 + var / mu ** 2); mln = np.log(mu) - 0.5 * s2; s = np.sqrt(s2)
    cdf = lambda x: lognorm.cdf(x, s=s, scale=np.exp(mln))
    k = np.arange(0, K + 1); lo = np.maximum(k - 0.5, 0.0)
    return cdf(k + 0.5) - cdf(lo)

def tv(p, q):
    n = max(len(p), len(q)); P, Q = np.zeros(n), np.zeros(n)
    P[:len(p)] = p; Q[:len(q)] = q
    return 0.5 * np.abs(P - Q).sum()

def kl(p, q, eps=1e-12):
    n = max(len(p), len(q)); P, Q = np.zeros(n), np.zeros(n)
    P[:len(p)] = p; Q[:len(q)] = q; mask = P > eps
    return float(np.sum(P[mask] * np.log(P[mask] / np.clip(Q[mask], eps, None))))

def central_interval(p, level=0.90):
    """Smallest central interval [lo,hi] with cdf mass >= level; return (lo,hi)."""
    cdf = np.cumsum(p); cdf /= cdf[-1]
    lo = np.searchsorted(cdf, (1 - level) / 2)
    hi = np.searchsorted(cdf, 1 - (1 - level) / 2)
    return lo, hi

def coverage_under_exact(pe, psurr, level=0.90):
    """Exact probability mass falling inside the surrogate's nominal-`level`
    central interval (calibration of the surrogate interval)."""
    lo, hi = central_interval(psurr, level)
    return pe[lo:hi + 1].sum()

def part1():
    print("=== Part 1: low-count transition-approximation error (p_s=0.6, a set to hit mu) ===")
    ps = 0.6
    # choose (n, a) so that mu = n*ps + a hits the target; keep a modest so survival matters
    targets = [1.0, 1.5, 2.0, 3.0, 5.0, 9.6]
    print(f"{'mu':>5} {'TVg':>7} {'TVln':>7} {'KLg':>7} {'KLln':>7} "
          f"{'P0ex':>7} {'P0g':>7} {'P0ln':>7} {'cov90g':>7} {'cov90ln':>8}")
    rows = []
    for mu in targets:
        a = max(0.5, 0.4 * mu)           # arrival share
        n = max(1, int(round((mu - a) / ps)))
        a = mu - n * ps                  # exact arrival so mean == mu
        if a < 0:
            n = int(np.floor(mu / ps)); a = mu - n * ps
        var = n * ps * (1 - ps) + a; sig = np.sqrt(max(var, 1e-9))
        K = int(mu + 12 * sig + 60)
        pe = exact_pmf(n, ps, a, K)
        pg = gaussian_discrete(mu, sig, K)
        pl = lognormal_discrete(mu, var, K)
        cg = coverage_under_exact(pe, pg); cl = coverage_under_exact(pe, pl)
        rows.append((mu, tv(pe, pg), tv(pe, pl), kl(pe, pg), kl(pe, pl),
                     pe[0], pg[0], pl[0], cg, cl))
        print(f"{mu:5.1f} {tv(pe,pg):7.3f} {tv(pe,pl):7.3f} {kl(pe,pg):7.3f} "
              f"{kl(pe,pl):7.3f} {pe[0]:7.3f} {pg[0]:7.3f} {pl[0]:7.3f} "
              f"{cg:7.3f} {cl:8.3f}")
    return rows

# ----------------------------------------------------------------------
# Part 2: latent-field recovery, exact vs Gaussian-transition particle filter
# ----------------------------------------------------------------------
def build_grid(side):
    coords = np.array([(i, j) for i in range(side) for j in range(side)], float)
    m = len(coords)
    W = np.zeros((m, m))
    for i in range(m):
        d = np.abs(coords - coords[i]).max(1)
        nb = (d == 1)
        W[i, nb] = 0.5 / max(nb.sum(), 1)      # normalized queen weights
    return coords, W

def nb_logpmf(y, mu, r):
    mu = np.maximum(mu, 1e-9)
    return (lgamma_vec(y + r) - lgamma_vec(r) - lgamma_vec(y + 1.0)
            + r * np.log(r / (r + mu)) + y * np.log(mu / (r + mu)))

_lg = np.vectorize(lgamma)
def lgamma_vec(x):
    return _lg(x)

def simulate_truth(side, T, phi, K, psi, eta, rho, r, seed):
    coords, W = build_grid(side); m = len(coords)
    rr = np.random.default_rng(seed)
    N = np.zeros((T, m), int)
    N[0] = rr.poisson(max(psi / (1 - phi), 1.0), m)
    for t in range(1, T):
        ps = phi * np.exp(-N[t - 1] / K)
        g = N[t - 1] / (N[t - 1] + eta)
        Lam = W @ (g * N[t - 1])
        N[t] = rr.binomial(N[t - 1], ps) + rr.poisson(Lam) + rr.poisson(psi, m)
    mu_obs = rho * N
    Y = rr.poisson(rr.gamma(r, mu_obs / r + 1e-12))       # NB(mean=rho N, disp r)
    return N, Y, W

def run_pf(Y, W, phi, K, psi, eta, rho, r, mode, P=4000, seed=0):
    """Bootstrap particle filter over the integer field. mode in {'exact','gauss'}."""
    T, m = Y.shape
    rr = np.random.default_rng(seed)
    # init particles from a diffuse Poisson around observed level
    init_mean = np.maximum(Y[0] / max(rho, 1e-3), 0.5)
    part = rr.poisson(np.broadcast_to(init_mean, (P, m))).astype(float)
    means = np.zeros((T, m)); lo = np.zeros((T, m)); hi = np.zeros((T, m))
    def record(t):
        srt = np.sort(part, axis=0)
        means[t] = part.mean(0)
        lo[t] = srt[int(0.05 * P)]; hi[t] = srt[min(int(0.95 * P), P - 1)]
    # weight t=0
    lw = nb_logpmf(Y[0][None, :], rho * part, r).sum(1)
    w = np.exp(lw - lw.max()); w /= w.sum()
    idx = systematic_resample(w, rr); part = part[idx]
    record(0)
    for t in range(1, T):
        ps = phi * np.exp(-part / K)
        g = part / (part + eta)
        Lam = (g * part) @ W.T                          # (P,m)
        if mode == 'exact':
            S = rr.binomial(part.astype(int), ps)
            I = rr.poisson(Lam)
            L = rr.poisson(psi, size=(P, m))
            part = (S + I + L).astype(float)
        elif mode == 'gauss':  # gauss moment-matched transition, rounded & clipped at 0
            mu_c = part * ps + Lam + psi
            var_c = part * ps * (1 - ps) + Lam + psi
            draw = rr.normal(mu_c, np.sqrt(np.maximum(var_c, 1e-9)))
            part = np.clip(np.round(draw), 0, None)
        elif mode == 'lnorm':  # moment-matched log-normal relaxation (continuous, positive)
            mu_c = part * ps + Lam + psi
            var_c = part * ps * (1 - ps) + Lam + psi
            mu_c = np.maximum(mu_c, 1e-3)
            s2 = np.log(1.0 + var_c / mu_c ** 2)
            draw = np.exp(rr.normal(np.log(mu_c) - 0.5 * s2, np.sqrt(s2)))
            part = draw
        lw = nb_logpmf(Y[t][None, :], rho * part, r).sum(1)
        w = np.exp(lw - lw.max()); w /= w.sum()
        idx = systematic_resample(w, rr); part = part[idx]
        record(t)
    return means, lo, hi

def systematic_resample(w, rr):
    P = len(w); positions = (rr.random() + np.arange(P)) / P
    cum = np.cumsum(w); cum[-1] = 1.0
    return np.searchsorted(cum, positions)

def part2():
    print("\n=== Part 2: latent-field recovery, exact vs Gaussian vs log-normal transition PF ===")
    regimes = {
        'low':  dict(phi=0.30, K=6.0,  psi=0.8, eta=3.0, rho=0.6, r=6.0),
        'app':  dict(phi=0.45, K=40.0, psi=4.0, eta=8.0, rho=0.6, r=6.0),
    }
    side, T = 4, 8
    print(f"{'regime':8} {'meanN':>6} {'filter':>8} {'RMSE':>7} {'cov90':>7}")
    results = {}
    for name, pr in regimes.items():
        acc = {m: dict(rmse=[], cov=[]) for m in ('exact', 'gauss', 'lnorm')}
        mns = []
        for rep in range(5):
            N, Y, W = simulate_truth(side, T, seed=300 + rep, **pr)
            mns.append(N[1:].mean()); Nt = N[1:]
            for mode in ('exact', 'gauss', 'lnorm'):
                mm, lo, hi = run_pf(Y, W, mode=mode, seed=10 + rep, **pr)
                acc[mode]['rmse'].append(np.sqrt(((mm[1:] - Nt) ** 2).mean()))
                acc[mode]['cov'].append(((Nt >= lo[1:]) & (Nt <= hi[1:])).mean())
        mn = np.mean(mns); results[name] = dict(meanN=mn)
        for mode in ('exact', 'gauss', 'lnorm'):
            rmse = np.mean(acc[mode]['rmse']); cov = np.mean(acc[mode]['cov'])
            results[name][mode] = dict(rmse=rmse, cov=cov)
            print(f"{name:8} {mn:6.2f} {mode:>8} {rmse:7.3f} {cov:7.3f}")
    return results

if __name__ == "__main__":
    part1()
    part2()
