# Data

## Processed arrays (included)
`mos_arrays.npz` — NEON UNDE 2024 *Coquillettidia perturbans* dual-channel
arrays built by `code/build_mos_arrays.py`: species-specific overnight counts
(`YC`), species-specific daytime detections (`YP`), effort, land-cover/elevation
covariates, plot coordinates, and the 2.5 km dispersal adjacency (`W`).

## Raw source (not redistributed here)
NEON Mosquitoes sampled from CO2 traps, product **DP1.10043.001**, RELEASE-2024,
site UNDE, 2024 season. Openly available from the NEON data portal
(https://data.neonscience.org). Place the stacked tables
`mos_trapping_stacked.csv`, `mos_sorting_stacked.csv`,
`mos_expertTaxonomistIDProcessed_stacked.csv` alongside `build_mos_arrays.py`
and rerun it to regenerate `mos_arrays.npz`.

NEON data are provided under a CC0 public-domain dedication.
