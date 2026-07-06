"""
Joint-estimation robustness study for SC-ABM-NKD (addresses reviewer points
3.2/3.3: whether the abundance scale is data-driven or prior-driven under
count-only surveys).

Supplement study C is the CONSERVATIVE upper bound on the detection-abundance
confound: it conditions on the TRUE latent field and scales carrying capacity
with the abundance factor, so count-only data see only the product
expit(mu_rho)*c*N and the scale is a pure likelihood ridge broken only by the
prior. This study restores the dynamics. It jointly estimates the detection
intercept mu_rho AND the carrying capacity K (hence the absolute abundance
scale) from a count-only series, with the remaining process parameters
(phi, psi, r) anchored at their values --- the regime in which S1/S3 show the
process parameters are recovered. The density-dependent approach to
carrying capacity gives the count series a shape that pins the absolute level,
so mu_rho becomes identified rather than prior-driven.

Inference is an EXACT hidden-Markov forward filter on a truncated integer
state space (no Monte Carlo): the marginal likelihood p(Y_{1:T} | mu_rho, K)
is computed exactly, so the posterior is exact up to the state truncation.
Self-contained (numpy + scipy).
"""
import numpy as np
from math import lgamma
from scipy.stats import binom, poisson

NMAX = 200
STATES = np.arange(NMAX + 1)
_lg = np.vectorize(lgamma)

# ---- true data-generating process (single density-regulated population) ----
PHI, KTRUE, PSI, RDISP = 0.60, 40.0, 6.0, 6.0
MU_TRUE = -0.5                      # detection intercept; rho = expit(-0.5) = 0.378
RHO_TRUE = 1 / (1 + np.exp(-MU_TRUE))
T = 16
N0 = 4

def simulate(seed):
    rr = np.random.default_rng(seed)
    N = np.zeros(T + 1, int); N[0] = N0
    for t in range(1, T + 1):
        ps = PHI * np.exp(-N[t - 1] / KTRUE)
        N[t] = rr.binomial(N[t - 1], ps) + rr.poisson(PSI)
    N = N[1:]
    mu_obs = RHO_TRUE * N
    Y = rr.poisson(rr.gamma(RDISP, mu_obs / RDISP + 1e-12))    # NB(mean=rho N, r)
    return N, Y

# ---- exact transition matrix for a given K ----
def transition_matrix(K):
    Tm = np.zeros((NMAX + 1, NMAX + 1))
    pois = poisson.pmf(STATES, PSI)
    for n in range(NMAX + 1):
        ps = PHI * np.exp(-n / K)
        b = binom.pmf(np.arange(n + 1), n, ps)
        conv = np.convolve(b, pois)[:NMAX + 1]
        row = np.zeros(NMAX + 1); row[:len(conv)] = conv
        s = row.sum()
        row[NMAX] += max(0.0, 1.0 - s)          # pile truncated tail on NMAX
        Tm[n] = row
    return Tm

def nb_obs_matrix(mu_rho):
    """Obs[t, n] = P(Y_t | N=n, mu_rho)."""
    rho = 1 / (1 + np.exp(-mu_rho))
    mean = rho * STATES + 1e-9
    logC = _lg(Yobs[:, None] + RDISP) - _lg(RDISP) - _lg(Yobs[:, None] + 1.0)
    logp = logC + RDISP * np.log(RDISP / (RDISP + mean))[None, :] \
              + Yobs[:, None] * np.log(mean / (RDISP + mean))[None, :]
    return np.exp(logp)

def forward_loglik(Tm, Obs, p0):
    """Exact HMM marginal log-likelihood."""
    alpha = p0 * Obs[0]
    ll = np.log(alpha.sum()); alpha /= alpha.sum()
    for t in range(1, T):
        alpha = (alpha @ Tm) * Obs[t]
        c = alpha.sum(); ll += np.log(c); alpha /= c
    return ll

def expected_meanN(Tm, p0):
    """Exact prior-predictive time-averaged E[N] under dynamics with this K."""
    p = p0.copy(); tot = 0.0
    for t in range(T):
        p = p @ Tm
        tot += (STATES * p).sum()
    return tot / T

# ---------------------------------------------------------------------------
def loggamma_logpdf(x, a, b):
    return a * np.log(b) - lgamma(a) + (a - 1) * np.log(x) - b * x

def run(K_prior='anchored', nrep=6):
    global Yobs
    mus = np.linspace(-3.5, 3.5, 141)
    Ks = np.linspace(12, 130, 60)
    p0 = np.zeros(NMAX + 1); p0[N0] = 1.0

    Tmats = [transition_matrix(K) for K in Ks]
    meanN_of_K = np.array([expected_meanN(Tm, p0) for Tm in Tmats])

    # log-likelihood surface averaged over replicate datasets (reduces DGP noise)
    LL = np.zeros((len(mus), len(Ks)))
    for rep in range(nrep):
        Ntrue, Yobs = simulate(1000 + rep)
        Obs_by_mu = [nb_obs_matrix(mu) for mu in mus]
        for j, Tm in enumerate(Tmats):
            for i, Obs in enumerate(Obs_by_mu):
                LL[i, j] += forward_loglik(Tm, Obs, p0)
    LL /= nrep

    # K prior
    if K_prior == 'anchored':
        lKpr = -0.5 * ((np.log(Ks) - np.log(KTRUE)) / 0.35) ** 2      # informative
    else:  # diffuse main-analysis-like prior on log K
        lKpr = -0.5 * ((np.log(Ks) - np.log(KTRUE)) / 1.5) ** 2

    print(f"\n=== Joint estimation of (mu_rho, K), count-only, K-prior={K_prior} ===")
    print(f"truth: mu_rho={MU_TRUE}, rho={RHO_TRUE:.3f}, K={KTRUE}, mean N_true~"
          f"{meanN_of_K[np.argmin(abs(Ks-KTRUE))]:.1f}, prior sd(mu_rho)=see rows")
    print(f"{'mu_rho prior':22s}{'post mean':>10}{'post sd':>9}{'rho_bar':>9}"
          f"{'Khat':>8}{'Nhat':>8}")
    table = []
    for name, (m_pr, s_pr) in {
        "N(0,1) [main]":     (0.0, 1.0),
        "N(-0.5,1) [truth]": (-0.5, 1.0),
        "N(0,0.5^2) [tight]":(0.0, 0.5),
        "N(-1.5,1)":         (-1.5, 1.0),
        "N(0.5,1)":          (0.5, 1.0),
    }.items():
        lmupr = -0.5 * ((mus - m_pr) / s_pr) ** 2
        logpost = LL + lmupr[:, None] + lKpr[None, :]
        logpost -= logpost.max()
        P = np.exp(logpost); P /= P.sum()
        pmu = P.sum(1); pK = P.sum(0)
        mu_m = (mus * pmu).sum()
        mu_sd = np.sqrt(((mus - mu_m) ** 2 * pmu).sum())
        rho_m = ((1 / (1 + np.exp(-mus))) * pmu).sum()
        K_m = (Ks * pK).sum()
        N_m = (meanN_of_K * pK).sum()
        table.append((name, s_pr, mu_m, mu_sd, rho_m, K_m, N_m))
        print(f"{name:22s}{mu_m:10.2f}{mu_sd:9.2f}{rho_m:9.3f}{K_m:8.1f}{N_m:8.2f}")
    return table

if __name__ == "__main__":
    run('anchored')
    run('diffuse')
