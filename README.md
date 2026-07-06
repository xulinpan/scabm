# SC-ABM: Multi-Source Data Fusion and Observation-Layer Identifiability

Reproducibility repository for the manuscript

> **Multi-Source Data Fusion and Observation-Layer Identifiability for
> Count-Valued Spatio-Temporal Abundance Models**
> (submitted to *Journal of Agricultural, Biological and Environmental Statistics*)

A count-valued spatio-temporal agent-based / state-space model observed through
two conditionally independent channels — a negative-binomial count survey and a
false-positive-aware Bernoulli detection survey — with an analytic
Fisher-information theory of when the two streams jointly identify the
observation parameters, a value-of-information simulation, and a real
dual-channel NEON mosquito application.

## Layout

```
manuscript/   LaTeX source, figures, response-to-reviewers, changelog
code/         model, estimators, simulation studies, NEON pipeline
data/         processed NEON arrays (mos_arrays.npz)
results/      posterior summaries (neon_post_summary.json)
```

## Reproduce the NEON application

```bash
pip install -r requirements.txt
python code/build_mos_arrays.py        # build mos_arrays.npz from raw NEON tables
python code/mos_numpyro_fit.py         # NUTS fit -> mos_numpyro_posterior.npz + summary
python code/make_mos_figs.py           # regenerate mos_field.png, mos_ppc.png
```

Set `FREE_RHO=1` before `mos_numpyro_fit.py` to estimate the detection fraction
rho under an informative prior instead of fixing rho=1.

## Reproduce the simulation studies

```bash
python code/paper2_sim_studies.py           # coverage studies (checkpointed)
python code/paper2_sim_studies.py summarize  # coverage / width / bias table
```

## Key result

On the corrected species-specific NEON channel the detection intensity is
well identified (kappa = 0.125 [0.065, 0.222], split-R-hat 1.00), while the
absolute abundance scale is anchored by rho=1 and only weakly identified
(peak total abundance 32,600 [23,500, 56,500]).

## Data

NEON mosquito data (product DP1.10043.001) are openly available from the NEON
data portal; see `data/README.md`. The processed arrays used by the code are
included in `data/mos_arrays.npz`.

## License

Code released under the MIT License (see `LICENSE`). The manuscript text and
figures are © the authors.
