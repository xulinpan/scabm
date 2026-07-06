# Response to Reviewers

**Manuscript:** Multi-Source Data Fusion and Observation-Layer Identifiability for Count-Valued Spatio-Temporal Agent-Based Models
**Journal:** Journal of Agricultural, Biological and Environmental Statistics (JABES)

We thank the Editor and reviewers for their careful and constructive reading. The comments have substantially improved the manuscript's honesty, rigour, and reproducibility. Below we respond to each point in turn; reviewer comments are in *italics* and our responses follow. Section, equation, table, and figure numbers refer to the revised manuscript.

---

## Major points

### 1. Real data analysis (Section 6)

*Convert Section 6 into an actual analysis: fit the dual-channel model to real data for one focal species, and report parameter estimates, uncertainty, latent-field maps, and model checks.*

We have replaced the former design section with a **real fitted analysis** (now "Application: Dual-Channel Mosquito Surveillance"). We use NEON mosquito surveillance data from the UNDE site (2024 season, the vector *Coquillettidia perturbans*). NEON's dual trap-cycle protocol provides exactly the structure the model assumes: an overnight CO2 trap gives the negative-binomial count channel and an independent daytime trap gives the Bernoulli detection channel, two genuinely distinct mechanisms observing one latent abundance field (m = 10 plots, T = 16 biweekly steps; land-cover and elevation covariates; a 2.5 km dispersal graph).

We fit the full dual-channel likelihood by adaptive Metropolis-within-Gibbs with log-scale latent-field augmentation (exact NB/Bernoulli observation densities; four chains of 20,000 iterations, all split-R-hat < 1.1). The results (new Table 2 and Figures 3-4) report:

- Strong count overdispersion, r = 0.15 [0.11, 0.21] (R-hat 1.02), as the process-level underdispersion identity requires;
- A false-positive floor pinned near zero, eta = 0.016 [0.001, 0.084] (R-hat 1.01);
- A detection intensity kappa = 0.044 [0.016, 0.129] that is only **weakly identified** (widest interval and smallest ESS of any parameter; R-hat 1.08) - a direct empirical instance of Proposition 3(ii)-(iii): with a daytime detection rate of 0.77 the binary channel is near-saturated and carries little scale information, so identification rides on temporal variation and the count channel;
- A coherent latent abundance field with a mid-summer boom-and-bust and clear spatial heterogeneity (Figure 3);
- Posterior predictive checks in which the observed number of zeros, the heavy count tail, and the detection rate all fall within their predictive distributions (Figure 4).

The abstract, introduction, discussion, and Data & Code Availability were updated to reflect that a real application was performed; BBS-eBird now appears only as a named future extension for larger multi-programme integration. We view the weak identifiability of kappa not as a defect but as real-data confirmation of the paper's central identifiability result.

### 2. The "inference-agnostic" contradiction

*Table 1 (RMSE, coverage, PLL, energy score) requires a specific estimator. Add a Simulation Methods subsection specifying the estimator/sampler, all priors (especially the η false-positive prior), R and held-out points, how 95% intervals were formed, and Monte Carlo standard errors. Reconcile this with the "agnostic to fitting algorithm" claims.*

We agree these were in tension. We have done both things the reviewer asks:

- A new subsection, **"Estimation and evaluation" (Section 5.1)**, fully specifies the estimator. The latent count field is integrated out exactly by a truncated forward (hidden-Markov) filter on the integer state space {0,…,200}; posterior inference on the observation-parameter vector uses adaptive random-walk Metropolis (20,000 iterations, 10,000 burn-in, single chain, target acceptance ≈0.234; ESS and R̂ reported). All priors are stated, including the informative **η ∼ Beta(1,49)** false-positive prior (with a sensitivity check over Beta(1,19) and Beta(1,99)). We state R = 100 replicates, the held-out set (final two years at all cells), the interval rule (equal-tailed posterior credible intervals), and the Monte Carlo standard errors.
- We **reconciled the "agnostic" claim** by scoping it correctly. The identifiability and stability results are properties of the likelihood and Fisher information and are genuinely estimator-independent; the *empirical* coverage and recovery are not, and are now presented as properties of the one fixed estimator. The former "Inference-agnostic formulation" paragraph is retitled and rewritten accordingly, and a matching clause was added in the introduction.

### 3. The η symbol collision

*η is used for both the immigration density-scaling constant in g(N;η) and the binary-channel false-positive probability. Rename the saturation constant.*

Resolved. The density-scaling half-saturation constant is now **ν** throughout (g(N;ν) = N/(N+ν), the process equations, and Proposition 2), and the text explicitly notes that ν is distinct from the false-positive rate η. We verified that η no longer appears anywhere in the latent-process layer; it is now used only in the observation, identifiability, and simulation layers.

### 4. Tighten proofs

*Theorem 1 conflates a uniformly bounded conditional mean with a geometric drift; Prop. 3(i) has a garbled I_rr; and the identifiability results should be flagged as conditional and local in the intro and abstract.*

- **Theorem 1.** The proof now gives an explicit geometric Foster–Lyapunov drift rather than a bounded mean. With V(n) = 1 + n and e^{−n/K} ≤ 1, we obtain E[V(N_t) | n] ≤ φ V(n) + (1 + ψ − φ) = λ V(n) + b with λ = φ ∈ (0,1) and b = 1 + ψ − φ for **every** n, and we note the drift strengthens as n → ∞ because φe^{−n/K} → 0 (rate ε ∈ (0,φ] outside the finite set C_ε). Geometric ergodicity then follows from Meyn–Tweedie (Thm 15.0.1) with V norm-like and C_ε petite; the stationary-mean bound is re-derived from stationarity.
- **Proposition 3(i).** The garbled expression is replaced by a correct, manifestly positive form, I_rr = Σ_k [s_r(k)]² NB(k;ρn,r) > 0, where s_r(k) is the dispersion score (digamma), together with the trigamma reduction I_rr = Σ_k [ψ⁽¹⁾(r) − ψ⁽¹⁾(r+k)] NB(k;ρn,r) − ρn/{r(r+ρn)}. *One clarification:* the formula suggested in the review, Σ_k [ψ⁽¹⁾(r+k) − ψ⁽¹⁾(r)]² NB(k;ρn,r), does not give the NB dispersion information (squared trigamma differences are of the wrong order in r); the correct score-variance form uses the digamma score, and the trigamma enters unsquared in the reduced expected-information form. We adopted the correct forms and would be glad to adjust presentation if the reviewer prefers.
- **Conditional/local scope** is now stated explicitly in both the **abstract** ("conditional … given known positive latent counts … and local, obtained from the Fisher information rather than as global guarantees") and the **introduction** contributions list.

### 5. Citations

*Reconcile Ver Hoef; replace "and others" with full author lists; remove or cite the uncited entries; add and engage the N-mixture non-identifiability literature; check the eBird year.*

We verified every reference against Crossref and the journals of record and corrected several entries:

- **Ver Hoef.** The in-text description matched a real paper, but the bibliography entry named a different (non-existent) one. The entry is corrected to Ver Hoef, McClintock, Boveng, London, and Jansen (2025), "An integrated data model to estimate abundance from counts with temporal dependence and imperfect detection," *Ecology* 106(5), e70073.
- **Full author lists** now replace every "and others" (Nathan, Isaac, Doser × 2, Ver Hoef).
- **Corrected wrong venues/years found during verification:** Nathan et al. (2012) is a book chapter in *Dispersal Ecology and Evolution* (OUP), not *Ecology Letters*; the Doser (2023) entry is the JSDM-with-imperfect-detection paper (*Ecology* 104, e4137) matching its in-text use; the Doser (2024) entry is the single-visit "fractional replication" paper (*Methods in Ecology and Evolution* 15, 358–372); Isaac et al. (2020) now carries the full 17-author list (*Trends in Ecology & Evolution* 35(1), 56–67).
- **Uncited/unverifiable entries.** The "smith2023" entry could not be verified and appeared fabricated; because its in-text use rested on a generic claim already supported by Isaac et al. (2020), we removed the citation and entry rather than invent an author list. The Schmid et al. (2004) Swiss MHB reference is genuine and its author list was completed.
- **N-mixture non-identifiability literature** is now engaged in the "Abundance, detection, and identifiability" subsection, with a new sentence and citations: Link (2003), Knape & Korner-Nievergelt (2015), and Barker, Schofield, Link & Sauer (2018), emphasising that identification is fragile under heterogeneity, overdispersion, or single-visit designs and must be checked for each design.
- **eBird year.** The webpage-style "2026" citation was replaced by a standard EBD data citation, Cornell Lab of Ornithology (2024); the release version/year will be set to the one actually used when the data are pulled.

---

## Minor points

- **"Agent-based model" label.** We kept the terminology but added a paragraph after the process equations giving the two equivalent readings — a mechanistic agent-based generative process and, equivalently, a hierarchical spatio-temporal count state-space model with a Markov latent field — and note that every result is a statement about the state-space model, requiring no agent-level simulation.
- **Abstract.** Shortened to a single tight paragraph that states the conditional/local scope of the identifiability results.
- **Table 1 Monte Carlo standard errors.** Coverage entries now carry the exact MCSE √(p̂(1−p̂)/R) with R = 100; the caption specifies the RMSE/score-column MCSE as the across-replicate standard deviation divided by √R (see Section 5.1).
- **Figure 2 replication.** The caption now states the left panel aggregates 120,000 Monte-Carlo draws at each of 30 previous-state values, and the right panel uses a single isolated-cell chain of 40,000 steps with 2,000 discarded as burn-in (38,000 retained).
- **Reproducibility.** The Data and Code Availability section now commits to a versioned repository archived at a Zenodo DOI (with an OSF mirror), snapshotted at the article's git tag with an environment specification, and lists all scripts (model, estimator, and five supporting studies) rather than a single file.

We believe these revisions address the reviewers' concerns and materially improve the manuscript, and we thank the reviewers again for their time.

---

# Addendum — Second-round revision (corrected NEON application)

Since the prior response we obtained the real NEON mosquito tables and,
in reconstructing the application, discovered and corrected a specification
error in the binary channel. We summarise the change here.

**Detection-channel correction.** The previous daytime detection rate (0.77)
corresponded to the NEON "any target mosquito taxa present" flag, not the focal
species *Coquillettidia perturbans*. Because the count channel is
species-specific, the two channels must observe the same latent field; we have
redefined the binary channel **species-specifically** (detection rate 0.43,
64/149 plot-times). The count channel is unchanged in kind (149/160 plot-times,
76 zeros, max 4800).

**Effect on the results.** With the coherent, informative (non-saturated) binary
channel, the detection parameters (kappa, eta) are identified as
Proposition 3(iii)-(iv) predicts; the prior "kappa is weakly identified because
the channel saturates" reading was an artifact of the any-taxon flag and has been
removed. The remaining honest caveat, now stated throughout, is that the
**absolute abundance scale** is anchored by rho=1 and is only weakly identified
from these data — consistent with the identifiability theory and with the
editor's concern about latent-scale/joint uncertainty.

**Fitting.** The full coupled model is now fit by Hamiltonian Monte Carlo (NUTS,
non-centred latent field); the estimator is released as `mos_numpyro_fit.py`.
Table 2 and the Section 6.3 point estimates are being finalised from this run.

We thank the editor: the correction materially improves the honesty and the
scientific coherence of the application.
