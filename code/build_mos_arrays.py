"""
Reconstruct the NEON UNDE 2024 dual-channel mosquito arrays for
Coquillettidia perturbans from the stacked NEON tables (DP1.10043.001):

  count channel  Y^(C): overnight (night) CO2-trap counts, subsampling-corrected
  binary channel Y^(P): daytime-trap detection (0/1)

Joins: trapping(sampleID) -> sorting(sampleID,subsampleID,proportionIdentified)
       -> expertID(subsampleID, scientificName, individualCount).
Outputs mos_arrays.npz with YC, YP, effort, covariates, adjacency.
"""
import csv, numpy as np
from collections import defaultdict
from datetime import datetime

BASE = "/sessions/confident-clever-feynman/mnt/uploads"
SITE = "UNDE"; YEAR = "2024"; SP = "Coquillettidia perturbans"
START = datetime(2024, 4, 11); NBINS = 16; BIN = 14  # biweekly

def rd(f): return list(csv.DictReader(open(f"{BASE}/{f}")))
def cdate(s): return datetime.strptime(s[:10], "%Y-%m-%d")
def tbin(s):
    d = (cdate(s) - START).days
    return d // BIN if 0 <= d < NBINS * BIN else None

trap = [r for r in rd("mos_trapping_stacked.csv") if r['siteID']==SITE and r['collectDate'][:4]==YEAR]
sort = [r for r in rd("mos_sorting_stacked.csv")   if r['siteID']==SITE and r['collectDate'][:4]==YEAR]
exp  = [r for r in rd("mos_expertTaxonomistIDProcessed_stacked.csv") if r['siteID']==SITE and r['collectDate'][:4]==YEAR]

plots = sorted(set(r['plotID'] for r in trap)); m = len(plots); pidx = {p:i for i,p in enumerate(plots)}

# subsampleID -> proportionIdentified ; and subsampleID -> sampleID
prop = {}; sub2sample = {}
for r in sort:
    try: prop[r['subsampleID']] = float(r['proportionIdentified'])
    except: pass
    sub2sample[r['subsampleID']] = r['sampleID']

# sampleID -> corrected C.perturbans count (sum over its subsamples)
count_by_sample = defaultdict(float)
present_by_sample = defaultdict(bool)
for r in exp:
    ss = r['subsampleID']; samp = sub2sample.get(ss)
    if samp is None: continue
    if SP in r['scientificName']:
        n = float(r['individualCount']) if r['individualCount'] else 0.0
        p = prop.get(ss, 1.0) or 1.0
        count_by_sample[samp] += n / p
        if n > 0: present_by_sample[samp] = True

# fill channel arrays; night -> count, day -> binary
YC = np.full((m, NBINS), np.nan); EC = np.full((m, NBINS), np.nan)
YP = np.full((m, NBINS), np.nan); EP = np.full((m, NBINS), np.nan)
cov = {}
for r in trap:
    b = tbin(r['collectDate']);
    if b is None or r['plotID'] not in pidx: continue
    i = pidx[r['plotID']]; samp = r['sampleID']
    try: th = float(r['trapHours'])
    except: th = np.nan
    if th is not None and th == 0: continue
    nod = r['nightOrDay']
    if nod == 'night':
        c = count_by_sample.get(samp, 0.0)
        YC[i, b] = round(c) if np.isnan(YC[i,b]) else YC[i,b] + round(c)
        EC[i, b] = th
    elif nod == 'day':
        det = 1 if present_by_sample.get(samp, False) else 0
        YP[i, b] = det if np.isnan(YP[i,b]) else max(YP[i,b], det)
        EP[i, b] = th
    cov.setdefault(i, dict(elev=r['elevation'], nlcd=r['nlcdClass'],
                           lat=r['decimalLatitude'], lon=r['decimalLongitude']))

# covariates
elev = np.array([float(cov[i]['elev']) if cov.get(i) and cov[i]['elev'] else np.nan for i in range(m)])
wetland = np.array([1.0 if cov.get(i) and 'Wetlands' in cov[i]['nlcd'] else 0.0 for i in range(m)])
forest  = np.array([1.0 if cov.get(i) and 'orest' in cov[i]['nlcd'] else 0.0 for i in range(m)])
lat = np.array([float(cov[i]['lat']) for i in range(m)]); lon = np.array([float(cov[i]['lon']) for i in range(m)])

# 2.5 km adjacency (haversine)
def hav(a,b,c,d):
    R=6371.0; import math
    p1,p2=math.radians(a),math.radians(c); dphi=math.radians(c-a); dl=math.radians(d-b)
    h=math.sin(dphi/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R*math.asin(min(1,h**0.5))
W=np.zeros((m,m))
for i in range(m):
    for j in range(m):
        if i!=j and hav(lat[i],lon[i],lat[j],lon[j])<=2.5: W[i,j]=1
deg = W.sum(1)

# ---- summary vs manuscript ----
nC = int(np.isfinite(YC).sum()); zC = int((YC==0).sum()); maxC = int(np.nanmax(YC))
nP = int(np.isfinite(YP).sum()); det = int(np.nansum(YP))
print(f"plots m={m}, T={NBINS}")
print(f"COUNT channel: observed {nC}/{m*NBINS} plot-times, {zC} zeros, max {maxC}")
print(f"BINARY channel: observed {nP} plot-times, {det} detections (rate {det/max(nP,1):.2f})")
print(f"mean degree (2.5km): {deg.mean():.2f}")
print(f"nlcd wetland plots: {int(wetland.sum())}, forest plots: {int(forest.sum())}")

np.savez("/sessions/confident-clever-feynman/mnt/papers/mos_arrays.npz",
         YC=YC, YP=YP, EC=EC, EP=EP, elev=elev, wetland=wetland, forest=forest,
         lat=lat, lon=lon, W=W, plots=np.array(plots))
print("saved mos_arrays.npz")
