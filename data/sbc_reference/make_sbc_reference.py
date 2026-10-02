"""Golden reference data for the spectral bin classifier (SBC).

Runs the authoritative Python SBC by D. Tripp (2023-08-31) and writes the
results that SHARPlib's SBC is tested against:

    python3 make_sbc_reference.py [--reference-dir DIR] [--output-dir DIR]

--reference-dir defaults to /workspace/papers/sbc_python and --output-dir to
this script's directory. The reference files must match these SHA-256 sums:

    sbc_alg_2023Aug31.py  ac920c643038b215772c8b2cf620466e0a645025f98e1532cc5db8a69ee1daa6
    run_sbc.py            28c2497a9d80d8edce67d594abd103f14414b6e2246bdcdf8abfdf2a6e84929d
    sample_data.csv       12e7b756a138a6fb74a65cde1b938e7f4dd06cafd5f7e1352d8c235c7245fe34

Provenance
----------
- The reference files are read, never written. classify is loaded from
  sbc_alg_2023Aug31.py after the source substitutions listed under Patches,
  applied to an in-memory copy. Only the function definitions wetbulb_calc,
  ptop_slice, and pre_class are loaded from run_sbc.py, unpatched; its top
  level reads a hard-coded CSV path and is never run.
- The driver is run_sbc.py's: wetbulb_calc, ptop_slice, pre_class, classify,
  on top-down arrays with heights AGL (surface height 0 m).
- Every input is float32: level data, DSD diameters and concentrations, the
  ice nucleation temperature, and the rime factor. The reference runs on
  exactly those values promoted to float64; temp_nuc is 267.15f promoted
  (267.1499939 K), not run_sbc.py's literal 267.15. RH in percent is
  100 * relh, computed in float64.
- Wet-bulb: wetbulb_calc runs in float64 on the promoted float32 T, Td, and p.
  Its result is rounded to float32 and stored, and ptop_slice, pre_class, and
  classify run on the stored float32 wet-bulb promoted to float64. SHARPlib is
  given the stored wet-bulb, so both sides see identical inputs.
- deld is passed to classify as run_sbc.py passes it (0.7, 0.1, or 0.6, per
  DSD). It enters only mass_psd = deld * N * velocity_ratio, and the ratio is
  always 1, so it cancels in every mass ratio. The generator reruns every core
  case with deld = 1 and fails unless category, supercooled-liquid height,
  crossings, branches, flags, and the profile are identical and the liquid
  fraction agrees within 1e-12.

Patches to classify (source substitution; each anchor must match exactly)
-------------------------------------------------------------------------
1. Return rainw, raini, and crossings alongside the original outputs, so the
   liquid fraction is the unrounded rainw / (raini + rainw). The rounded
   percentage the reference returns is checked against it.
2. Branch tags: each level appends the letter of the branch it takes
   (A :138, B :164, C :186, D :217, E :243, F :261, G :402).
3. Tice switch fired: set when the assignment at :407 runs.
4. Shared refreeze level moved: set when :414 changes refrz_lvl while another
   bin already has refrz_lvl_flag set.
5. Initial rime factor: :9 reads a keyword argument (default 1.0, as in the
   reference) instead of the literal 1.0, for the rime-5 variants.
6. near_discontinuity instrumentation: records the fw of each melting update
   before its clamp (:316) and of each refreeze update (:438), psd_ri and
   psd_rw of each melting (:350) and refreezing (:459) bin, and the snow
   volume at both of its tests (:330 and :332) with the threshold in use.

Files (rows sorted by case_id, then level, then bin)
----------------------------------------------------
Levels run from the surface up: level 0 is the surface (2 m) level.
MISSING is -9999.

levels.parquet
    case_id int32; level int32; pressure (Pa), height (m AGL), temperature (K),
    dewpoint (K), relh (fraction), wetbulb (K, from wetbulb_calc): float32.
cases.parquet
    case_id int32
    group string: named, sample, or corpus
    description string
    dsd_name string: a dsds.parquet name
    rime_factor float32; ice_nucleation_temperature float32 (K)
    cloud_top_height float32 (m AGL; MISSING if no cloud)
    cloud_top_level int32 (from the surface; MISSING if no cloud)
    precip_type int32: RA 1, SN 2, RASN 3, FZRA 4, PL 5, FZRAPL 6, RAPL 7;
        MISSING if no cloud
    liquid_fraction float64: unrounded rainw / (raini + rainw) for core cases;
        SN 0, FZRA 1, RA 1 from the pre-classifier; MISSING if no cloud
    supercooled_liquid_height float32 (m AGL): the reference slw_hgt for core
        cases (MISSING where it is never assigned); SN MISSING, FZRA 0, RA
        MISSING from the pre-classifier; MISSING if no cloud
    stage string: no_cloud, preclassifier, or core
    crossings int32: the reference's 0 C crossing count; MISSING unless core
    branches string: the branch letters used, in order A-G (e.g. "ABFG");
        empty unless core
    tnuc_switched, refreeze_level_moved, near_threshold, near_discontinuity:
        bool (see Flags)
    disc_bin int32, disc_level int32 (from the surface), disc_scope string
        (bin or column): the first near_discontinuity hit; -1, -1, and empty
        when unflagged
profiles.parquet
    case_id int32; level int32; bin int32; liquid_fraction float32. Every
    level and bin of every case: the reference water_fraction mapped to
    bottom-up levels, MISSING above the cloud top and for every level of a
    case that is not core.
dsds.parquet
    dsd_name string; bin int32; diameter float32 (mm); concentration float32
    (number per bin; only ratios matter).

DSDs (dsd_name), from run_sbc.py's construction, stored as float32
------------------------------------------------------------------
python_default   deld 0.7: D = np.arange(0.05, 1.85 + deld, deld), N from
                 np.interp(D, diameter_orig, psd_orig). 4 bins: 0.05, 0.75,
                 1.45, 2.15 mm; 55.1843, 146.647, 11.6891, 3.60886.
python_deld_0.1  deld 0.1, the same construction. 19 bins, 0.05-1.85 mm.
cpp_2.0.3        deld 0.6: the MRMS C++ 2.0.3 literals, 0.05, 0.65, 1.25,
                 1.85 mm; 55.1843, 206.606, 25.4924, 3.60886.

Case groups and case_id ranges
------------------------------
named (0-999)      Idealized profiles of at most 15 levels. Pressure is
                   1000 hPa * exp(-z / 8 km), rounded to whole hPa. T, Td,
                   and relh are as listed in NAMED_CASES. A coverage table
                   proves they reach every category, pre-classifier path,
                   branch, switch, cloud-top rule, DSD, rime 5, and Tice -10 C.
sample (1000)      sample_data.csv as it is (relh = rh / 100), python_default.
corpus (10000+)    Perturbations of the sample, seeded with CORPUS_SEED. Each
                   family adds to T and Td a uniform shift and a Gaussian nose
                   a * exp(-((z - zc) / w)^2), with its own ranges (see
                   CORPUS_FAMILIES). With probability 0.5 a dry layer from zd
                   to zd + depth has Td = T - U(11, 25) K and relh U(0.05,
                   0.38). Corpus cases use python_default with probability
                   0.6, cpp_2.0.3 0.2, python_deld_0.1 0.2; rime 5 with
                   probability 0.15; Tice 263.15 K with probability 0.15.
                   Candidates whose cloud top is the surface level (a
                   one-level column) are skipped.

No exact-threshold inputs
-------------------------
A float32 value exactly at a threshold compares differently in the float64
reference than in SHARPlib's float comparison. So no case may contain:
- a wet-bulb exactly 273.15f, 273.65f, 263.15f, or the case's Tice;
- a relh exactly 0.60f, 0.40f, or 0.80f;
- two levels of the column tied, in float or in float64, for nearest to
  3000 m AGL.
The generator nudges such a value (T and Td by 0.01 K, relh by 0.0001, the
upper tied height by 1 m), notes it in the description, and fails if any
remain.

Flags
-----
near_threshold: a core liquid fraction within 1e-4 of 0.15, 0.60, or 0.85.
near_discontinuity: at some level, some bin's deciding value lies within 1e-5
(relative) of the threshold of one of the reference's hard switches, on
either side:
- fw against 1: the melting update before its clamp (:316-318) and the
  refreeze update (:438). fw(K-1) == 1 decides the F shortcut (:279), the
  Tice switch (:406), and the G refreeze test (:412).
- the per-bin ice and liquid mass ratios psd_ri / (psd_ri + psd_rw) and
  psd_rw / (psd_ri + psd_rw) against 0.15 (:357-364, :380-387, :466-487).
- the snow volume against 1.81e-5 * rime^3.26 at both tests (:330-335).
disc_bin and disc_level give the first hit, top-down and then by bin.
disc_scope is column if a G level lies at or below disc_level, otherwise bin.
Named cases carry neither flag: a flagged named profile has T and Td shifted
by 0.01 K steps until it is clear, or the generator fails.

Comparison rules
----------------
1. Consistency, always: the returned category equals the reference decision
   tree applied to the returned liquid fraction, with the case's crossings and
   surface class (surface Tw > 273.15 K). Pre-classified and no-cloud cases
   match their category and floats exactly.
2. Category: exactly precip_type, except that a corpus case flagged
   near_threshold may return the category on the other side of that
   threshold. Rule 1 still holds.
3. Surface values: liquid fraction within 1e-4; supercooled-liquid height
   and cloud top on the same level.
4. Profiles: within 1e-3 for every bin and level, except in a
   near_discontinuity case from disc_level downward: only disc_bin is exempt
   for disc_scope bin, and every bin for disc_scope column.
5. A near_discontinuity corpus case that breaks rules 2-3, or rule 4 beyond
   its exemption, is listed, not failed. Such cases stay at most 1 % of the
   corpus. Named cases carry no flags, so every rule applies to them in full.

The profile tolerance is ten times the liquid-fraction one. Two terms of the
reference cancel and magnify float32 rounding. One is the melting heat near
saturation at 0 C (:308). The other is the refreezing numerator near -31 C,
where unknwn_xsi changes sign (:427-428). At such a level a float32 port
differs by up to about 1e-3. The reference's own profile moves more than that
when T or Td at the level changes by one float32 ulp.
"""

import argparse
import ast
import hashlib
import os
import time
import warnings
from dataclasses import dataclass, replace

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

F32 = np.float32
F64 = np.float64
MISSING = -9999

REFERENCE_SHA256 = {
    "sbc_alg_2023Aug31.py": "ac920c643038b215772c8b2cf620466e0a645025f98e1532cc5db8a69ee1daa6",
    "run_sbc.py": "28c2497a9d80d8edce67d594abd103f14414b6e2246bdcdf8abfdf2a6e84929d",
    "sample_data.csv": "12e7b756a138a6fb74a65cde1b938e7f4dd06cafd5f7e1352d8c235c7245fe34",
}

PRECIP_TYPE = {"RA": 1, "SN": 2, "RASN": 3, "FZRA": 4, "PL": 5, "FZRAPL": 6, "RAPL": 7}
PRECLASSIFIED = {"SN": (0.0, MISSING), "FZRA": (1.0, 0.0), "RA": (1.0, MISSING)}
TICE_DEFAULT = 267.15
TICE_ALT = 263.15
LF_THRESHOLDS = (0.15, 0.60, 0.85)
NEAR_THRESHOLD = 1e-4
NEAR_DISC_REL = 1e-5
CORPUS_CAP = 0.01
DELD_CHECK = 1.0

CORPUS_SEED = 20261002
# name, count, shift (K), nose amplitude (K), nose height (m), nose width (m)
CORPUS_FAMILIES = (
    ("aloft", 800, (-3.0, 3.0), (-4.0, 4.0), (800.0, 2500.0), (150.0, 1200.0)),
    ("surface", 1000, (-4.0, 1.0), (2.0, 9.0), (0.0, 300.0), (100.0, 600.0)),
    ("broad", 600, (-2.0, 6.0), (-6.0, 6.0), (0.0, 3000.0), (150.0, 1500.0)),
)

PSD_ORIG = [55.1843, 66.0695, 130.272, 154.556, 203.649, 171.814, 206.606, 146.647, 94.9404,
            79.4013, 61.0083, 35.6567, 25.4924, 16.2522, 11.6891, 7.49152, 3.60886]
DIAMETER_ORIG = [0.05 + 0.1 * i for i in range(len(PSD_ORIG))]


@dataclass(frozen=True)
class DSD:
    diameter: np.ndarray
    concentration: np.ndarray
    deld: float


def _python_dsd(deld):
    diameter = np.arange(0.05, 1.85 + deld, deld)
    return DSD(diameter.astype(F32), np.interp(diameter, DIAMETER_ORIG, PSD_ORIG).astype(F32), deld)


DSDS = {
    "python_default": _python_dsd(0.7),
    "python_deld_0.1": _python_dsd(0.1),
    "cpp_2.0.3": DSD(np.array([0.05, 0.65, 1.25, 1.85], F32),
                     np.array([55.1843, 206.606, 25.4924, 3.60886], F32), 0.6),
}


# ---------------------------------------------------------------- reference

def _anchor(src, old, new):
    if src.count(old) != 1:
        raise RuntimeError(f"patch anchor matched {src.count(old)} times: {old!r}")
    return src.replace(old, new)


def _after(src, old, line):
    return _anchor(src, old, old + line + "\n")


def _patch_classify(src):
    src = _anchor(src, "deld,particle_count,temp_nuc):\n",
                  "deld,particle_count,temp_nuc,rime_factor_init=1.0):\n"
                  "  _inst={'branch':[],'tnuc_switched':False,'refrz_moved':False,"
                  "'fw':[],'ratio':[],'vs':[]}\n")
    src = _anchor(src, "  rime_factor=1.0 #", "  rime_factor=rime_factor_init #")
    branches = (
        ("A", "     if K==0 and wbTemp1d[0]<=temp_nuc:\n", 7),
        ("B", "     elif 0<K<cross_index[0] and wbTemp1d[0]<=temp_nuc and wbTemp1d[K]<kelvin_frz:\n", 7),
        ("C", "     elif K==0 and wbTemp1d[0]>temp_nuc:\n", 7),
        ("D", "     elif 0<K<cross_index[0] and temp_nuc<wbTemp1d[0]<=kelvin_frz and "
              "wbTemp1d[K]>temp_nuc and ice_nuc_tops[0]>=cross_index[0]:\n", 9),
        ("E", "     elif 0<K<cross_index[0] and wbTemp1d[0]>kelvin_frz and wbTemp1d[K]>temp_nuc:\n", 7),
        ("F", "     elif wbTemp1d[K]>=kelvin_frz:\n", 7),
        ("G", "     else:\n       \n       refrz_detected_flag=True\n", 7),
    )
    for tag, cond, indent in branches:
        src = _after(src, cond, " " * indent + f"_inst['branch'].append('{tag}')")
    src = _after(src, "            temp_nuc=tnuc_alt#Change the ice nucleation temperature\n",
                 "            _inst['tnuc_switched']=True")
    src = _anchor(src, "                refrz_lvl=K-1\n",
                  "                if refrz_lvl!=K-1 and refrz_lvl_flag.max()==1: _inst['refrz_moved']=True\n"
                  "                refrz_lvl=K-1\n")
    src = _after(src, "*heat*layer_depth[K]\n",
                 "              _inst['fw'].append((K,J,water_fraction[K,J]))")
    src = _after(src, "                water_fraction[K,J]=mass_water[K,J]/mass_hydro[0,J]\n",
                 "                _inst['fw'].append((K,J,water_fraction[K,J]))")
    src = _after(src, "\n           psd_ri=3600./0.918*mass_total_ice\n",
                 "           _inst['ratio'].append((K,J,psd_ri,psd_rw))")
    src = _after(src, "\n             psd_ri=3600./0.918*mass_total_ice\n",
                 "             _inst['ratio'].append((K,J,psd_ri,psd_rw))")
    record_vs = "           _inst['vs'].append((K,J,volume_snow[K,J],1.81e-5*rime_factor**3.26))"
    src = _after(src, "           volume_snow[K,J]=343.*(rhs/rime_factor)**1.443\n", record_vs)
    src = _anchor(src, "           if volume_snow[K,J]>=1.81e-5*rime_factor**3.26: \n",
                  record_vs + "\n           if volume_snow[K,J]>=1.81e-5*rime_factor**3.26: \n")
    src = _anchor(src, "  return ptype,lwf,psd_ptype,water_fraction,slw_hgt",
                  "  return ptype,lwf,psd_ptype,water_fraction,slw_hgt,rainw,raini,crossings,_inst")
    return src


@dataclass(frozen=True)
class Reference:
    classify: object
    wetbulb_calc: object
    ptop_slice: object
    pre_class: object
    sample_csv: str


def load_reference(ref_dir):
    for name, digest in REFERENCE_SHA256.items():
        with open(os.path.join(ref_dir, name), "rb") as f:
            got = hashlib.sha256(f.read()).hexdigest()
        if got != digest:
            raise RuntimeError(f"{name}: SHA-256 {got} is not the 2023-08-31 reference")
    with open(os.path.join(ref_dir, "sbc_alg_2023Aug31.py")) as f:
        alg = _patch_classify(f.read())
    alg_ns = {}
    exec(compile(alg, "sbc_alg_2023Aug31.py (patched)", "exec"), alg_ns)
    with open(os.path.join(ref_dir, "run_sbc.py")) as f:
        driver = f.read()
    driver_ns = {"np": np}
    wanted = {"wetbulb_calc", "ptop_slice", "pre_class"}
    for node in ast.parse(driver).body:
        if isinstance(node, ast.FunctionDef) and node.name in wanted:
            exec(compile(ast.get_source_segment(driver, node), "run_sbc.py", "exec"), driver_ns)
            wanted.discard(node.name)
    if wanted:
        raise RuntimeError(f"run_sbc.py lacks {sorted(wanted)}")
    return Reference(alg_ns["classify"], driver_ns["wetbulb_calc"], driver_ns["ptop_slice"],
                     driver_ns["pre_class"], os.path.join(ref_dir, "sample_data.csv"))


# -------------------------------------------------------------------- cases

@dataclass(frozen=True)
class Case:
    group: str
    description: str
    dsd_name: str
    rime_factor: F32
    tice: F32
    pressure: np.ndarray
    height: np.ndarray
    temperature: np.ndarray
    dewpoint: np.ndarray
    relh: np.ndarray


def make_case(group, description, pressure, height, temperature, dewpoint, relh,
              dsd_name="python_default", rime_factor=1.0, tice=TICE_DEFAULT):
    arrays = [np.asarray(a, F64).astype(F32) for a in (pressure, height, temperature, dewpoint, relh)]
    return Case(group, description, dsd_name, F32(rime_factor), F32(tice), *arrays)


@dataclass
class Result:
    wetbulb: np.ndarray
    stage: str = "no_cloud"
    cloud_top_level: int = MISSING
    precip_type: str = ""
    liquid_fraction: float = float(MISSING)
    supercooled_liquid_height: float = float(MISSING)
    crossings: int = MISSING
    branch_by_k: tuple = ()
    water_fraction: np.ndarray = None
    psd_ptype: np.ndarray = None
    tnuc_switched: bool = False
    refreeze_level_moved: bool = False
    near_threshold: bool = False
    near_discontinuity: bool = False
    disc: tuple = (-1, -1, "")
    disc_switch: str = ""

    @property
    def branches(self):
        return "".join(sorted(set(self.branch_by_k)))


def wetbulb32(ref, case):
    tw = ref.wetbulb_calc(case.temperature.astype(F64), case.dewpoint.astype(F64),
                          case.pressure.astype(F64))
    return np.asarray(tw, F64).astype(F32)


def _first_discontinuity(inst):
    hits = []
    for k, j, fw in inst["fw"]:
        if abs(fw - 1.0) <= NEAR_DISC_REL:
            hits.append((k, j, "fw"))
    for k, j, ri, rw in inst["ratio"]:
        if ri + rw > 0:
            for ratio in (ri / (ri + rw), rw / (ri + rw)):
                if abs(ratio - 0.15) <= NEAR_DISC_REL * 0.15:
                    hits.append((k, j, "ratio"))
    for k, j, vs, threshold in inst["vs"]:
        if abs(vs - threshold) <= NEAR_DISC_REL * threshold:
            hits.append((k, j, "snow_volume"))
    return min(hits) if hits else None


def run_reference(ref, case, deld=None):
    """The run_sbc.py driver on one case. Arrays are bottom-up; the reference's are top-down."""
    tw = wetbulb32(ref, case)
    res = Result(wetbulb=tw)
    p, z, t, td, w = (a[::-1].astype(F64) for a in
                      (case.pressure, case.height, case.temperature, case.dewpoint, tw))
    rh = 100.0 * case.relh[::-1].astype(F64)
    p1, h1, t1, td1, w1, rh1, found = ref.ptop_slice(p, z, t, td, w, rh)
    if not found:
        return res
    res.cloud_top_level = len(p1) - 1
    tice = F64(case.tice)
    sbc_ptype = ref.pre_class(w1, h1, tice)
    if not pd.isna(sbc_ptype):
        res.stage = "preclassifier"
        res.precip_type = sbc_ptype
        res.liquid_fraction, res.supercooled_liquid_height = PRECLASSIFIED[sbc_ptype]
        return res
    dsd = DSDS[case.dsd_name]
    out = ref.classify(p1, h1, t1, td1, w1, rh1, len(dsd.diameter), dsd.diameter.astype(F64),
                       dsd.deld if deld is None else deld, dsd.concentration.astype(F64), tice,
                       rime_factor_init=F64(case.rime_factor))
    ptype, lwf, psd_ptype, water_fraction, slw_hgt, rainw, raini, crossings, inst = out
    lf = rainw / (raini + rainw)
    if round(lf * 100, 1) != lwf:
        raise RuntimeError(f"unrounded liquid fraction {lf} disagrees with the reference's {lwf} %")
    res.stage = "core"
    res.precip_type = ptype
    res.liquid_fraction = float(lf)
    res.supercooled_liquid_height = float(MISSING) if np.isnan(slw_hgt) else float(slw_hgt)
    res.crossings = int(crossings)
    res.branch_by_k = tuple(inst["branch"])
    res.water_fraction = water_fraction
    res.psd_ptype = psd_ptype
    res.tnuc_switched = inst["tnuc_switched"]
    res.refreeze_level_moved = inst["refrz_moved"]
    res.near_threshold = any(abs(lf - x) <= NEAR_THRESHOLD for x in LF_THRESHOLDS)
    hit = _first_discontinuity(inst)
    if hit is not None:
        k, j, switch = hit
        scope = "column" if "G" in res.branch_by_k[k:] else "bin"
        res.near_discontinuity = True
        res.disc = (j, res.cloud_top_level - k, scope)
        res.disc_switch = switch
    if len(res.branch_by_k) != res.cloud_top_level + 1:
        raise RuntimeError("branch tags do not cover the column")
    return res


# -------------------------------------------------- exact-threshold inputs

def threshold_violations(case, res):
    tw_thresholds = {F32(273.15), F32(273.65), F32(263.15), case.tice}
    relh_thresholds = {F32(0.60), F32(0.40), F32(0.80)}
    found = [("wetbulb", int(i)) for i in np.flatnonzero(np.isin(res.wetbulb, list(tw_thresholds)))]
    found += [("relh", int(i)) for i in np.flatnonzero(np.isin(case.relh, list(relh_thresholds)))]
    if res.cloud_top_level != MISSING:
        z = case.height[:res.cloud_top_level + 1]
        for dist in (np.abs(z - F32(3000.0)), np.abs(z.astype(F64) - 3000.0)):
            nearest = np.flatnonzero(dist == dist.min())
            if len(nearest) > 1:
                found.append(("height", int(nearest.max())))
    return sorted(set(found))


def nudge(case, violations):
    t, td, relh, z = (a.astype(F64) for a in (case.temperature, case.dewpoint, case.relh, case.height))
    for kind, i in violations:
        if kind == "wetbulb":
            t[i] += 0.01
            td[i] += 0.01
        elif kind == "relh":
            relh[i] += 0.0001
        else:
            z[i] += 1.0
    if np.any(np.diff(z) <= 0):
        raise RuntimeError(f"{case.description}: a height nudge broke monotonicity")
    note = ", ".join(f"{kind} at level {i}" for kind, i in violations)
    return replace(case, description=f"{case.description} [nudged: {note}]",
                   temperature=t.astype(F32), dewpoint=td.astype(F32), relh=relh.astype(F32),
                   height=z.astype(F32))


def settle(ref, case, max_rounds=10):
    """Runs the reference, nudging exact-threshold inputs until none remain."""
    for _ in range(max_rounds):
        res = run_reference(ref, case)
        violations = threshold_violations(case, res)
        if not violations:
            return case, res
        case = nudge(case, violations)
    raise RuntimeError(f"{case.description}: exact-threshold inputs remain")


# -------------------------------------------------------------- named cases

Z_COARSE = (0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000)
Z_FINE = (0, 250, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000)
Z_LOW = (0, 500, 1000, 1500, 2000, 2500, 3000)

SN_T = (0.2, -1, -2, -3, -4, -5, -7, -9, -11, -13, -17)
RASN_T = (0.45, -1, -2, -3, -4, -5, -7, -9, -11, -13, -17)
PL_T = (-5, -8, -7, -2, 0.4, -1, -3, -6.5, -9, -12, -17)
RAPL_T = (1, -1, -4, -8, -2, 0.4, -1, -3, -8, -12, -17)

# description, heights (m AGL), T, Td (None: Td = T), relh in % (None: 100), options.
# T and Td are in C unless kelvin=True.
NAMED_CASES = (
    ("SN from the pre-classifier: all subfreezing, cloud top and the 3 km layer colder than Tice",
     Z_COARSE, (-3, -4, -5, -7, -8, -9, -11, -13, -15, -17, -20), None, None, {}),
    ("FZRA from the pre-classifier: all subfreezing, cloud top warmer than Tice",
     Z_COARSE[:5], (-1, -2, -3, -4, -4.5), None, None, {}),
    ("FZRA from the pre-classifier: all subfreezing, nothing colder than Tice from 3 km down",
     Z_COARSE, (-1, -1.5, -2, -2.5, -3, -3.5, -4, -5.5, -8, -11, -14), None, None, {}),
    ("FZRA from the pre-classifier: cloud top warmer than Tice over a subfreezing surface",
     Z_COARSE[:6], (-2, -1, 2, 3, 1, -3), None, None, {}),
    ("RA from the pre-classifier: all warm",
     Z_COARSE[:5], (8, 6, 4, 3, 2), None, None, {}),
    ("RA from the core: frozen top, shallow warm layer at the surface",
     Z_FINE, (0.8, -1, -2, -3, -4, -5, -7, -9, -11, -13, -17), None, None, {}),
    ("SN from the core: frozen top, shallow warm layer at the surface",
     Z_FINE, SN_T, None, None, {}),
    ("RASN from the core: frozen top, shallow warm layer at the surface",
     Z_FINE, RASN_T, None, None, {}),
    ("PL from the core: partial melting, then refreezing that moves the shared refreeze level",
     Z_FINE, PL_T, None, None, {}),
    ("FZRAPL from the core: partial melting, then partial refreezing above Tice",
     Z_FINE, (-4, -5, -3, -1, 0.4, -1, -3, -6.5, -9, -12, -17), None, None, {}),
    ("FZRA from the core: full melting switches Tice to -10 C, so the -8 C layer does not refreeze",
     Z_FINE, (-3, -8, -8, -2, 2, 1, -1, -3, -8, -12, -17), None, None, {}),
    ("RAPL from the core: refrozen pellets melt again over a warm surface",
     Z_FINE, RAPL_T, None, None, {}),
    ("RA from the core: supercooled liquid cloud top over a warm surface",
     Z_LOW, (5, 4, 2, 1, -1, -2, -3), None, None, {}),
    ("RA from the core: warm cloud top above a subfreezing layer, warm surface",
     Z_LOW, (5, 3, -2, -3, 1, 2, 1.5), None, None, {}),
    ("RA from the core: a carried FZ class leaves the supercooled-liquid height at 2000 m",
     (0, 1000, 2000, 3000, 4000), (277, 268, 269, 280, 265), None, None, {"kelvin": True}),
    ("No cloud: too dry at every level",
     Z_COARSE[:5], (5, 3, 1, -1, -3), (-7, -9, -11, -13, -15), (30,) * 5, {}),
    ("Dry-layer restart: the cloud top moves below a dry layer at 3-3.5 km",
     Z_FINE, RASN_T, RASN_T[:7] + (-24, -26, -13, -17), (95,) * 7 + (25, 25, 95, 95), {}),
    ("RH >= 80 % fallback: no level passes the dewpoint test",
     Z_LOW, (-2, -4, -6.5, -8, -10.5, -13, -16), (-10, -12, -14.5, -16, -18.5, -21, -24),
     (50, 50, 50, 85, 85, 50, 50), {}),
    ("RASN from the core with the C++ 2.0.3 DSD",
     Z_FINE, RASN_T, None, None, {"dsd_name": "cpp_2.0.3"}),
    ("RAPL from the core with the deld 0.1 DSD (19 bins)",
     Z_FINE, RAPL_T, None, None, {"dsd_name": "python_deld_0.1"}),
    ("SN from the core with rime factor 5 (RASN with rime factor 1)",
     Z_FINE, RASN_T, None, None, {"rime_factor": 5.0}),
    ("FZRAPL from the core with Tice -10 C (PL with Tice -6 C)",
     Z_FINE, PL_T, None, None, {"tice": TICE_ALT}),
    ("SN from the core with the deld 0.1 DSD (19 bins)",
     Z_FINE, SN_T, None, None, {"dsd_name": "python_deld_0.1"}),
    ("RASN from the core with Tice -10 C (no refreezing, so Tice only decides the frozen top)",
     Z_FINE, RASN_T, None, None, {"tice": TICE_ALT}),
)


def named_pressure(z):
    return 100.0 * np.round(1000.0 * np.exp(-np.asarray(z, F64) / 8000.0))


def named_cases():
    cases = []
    for description, z, t, td, rh, opts in NAMED_CASES:
        opts = dict(opts)
        offset = 0.0 if opts.pop("kelvin", False) else 273.15
        t = np.asarray(t, F64) + offset
        td = t if td is None else np.asarray(td, F64) + offset
        relh = np.ones(len(z)) if rh is None else np.asarray(rh, F64) / 100.0
        if len(z) > 15:
            raise RuntimeError(f"{description}: more than 15 levels")
        cases.append(make_case("named", description, named_pressure(z), z, t, td, relh, **opts))
    return cases


def clear_named(ref, case):
    """Shifts T and Td of a flagged named profile by 0.01 K steps until no flag remains."""
    case, res = settle(ref, case)
    if not (res.near_threshold or res.near_discontinuity):
        return case, res
    for step in range(1, 41):
        offset = 0.01 * ((step + 1) // 2) * (1 if step % 2 else -1)
        shifted = replace(case, description=f"{case.description} [T, Td shifted {offset:+.2f} K]",
                          temperature=(case.temperature.astype(F64) + offset).astype(F32),
                          dewpoint=(case.dewpoint.astype(F64) + offset).astype(F32))
        shifted, res = settle(ref, shifted)
        if not (res.near_threshold or res.near_discontinuity):
            return shifted, res
    raise RuntimeError(f"{case.description}: cannot clear its flags")


# ------------------------------------------------------- sample and corpus

def sample_case(ref):
    df = pd.read_csv(ref.sample_csv).iloc[::-1]
    return make_case("sample", "HRRR 2022-02-02 23Z f00 at 34.94 N, 97.18 W (sample_data.csv)",
                     df["pres"], df["hgt_AGL"], df["temp"], df["dew"], df["rh"].to_numpy() / 100.0)


def corpus_candidates(base, family, count, rng):
    """Vectorized draw of `count` perturbed profiles of one family."""
    name, _, shift_range, amp_range, zc_range, width_range = family
    shift = rng.uniform(*shift_range, count)
    amp = rng.uniform(*amp_range, count)
    zc = rng.uniform(*zc_range, count)
    width = rng.uniform(*width_range, count)
    dry = rng.random(count) < 0.5
    dry_bottom = rng.uniform(-1000.0, 9000.0, count)
    dry_depth = rng.uniform(300.0, 15000.0, count)
    dry_depression = rng.uniform(11.0, 25.0, count)
    dry_relh = rng.uniform(0.05, 0.38, count)
    u = rng.random(count)
    dsd_names = np.where(u < 0.6, "python_default", np.where(u < 0.8, "cpp_2.0.3", "python_deld_0.1"))
    rime = np.where(rng.random(count) < 0.15, 5.0, 1.0)
    tice = np.where(rng.random(count) < 0.15, TICE_ALT, TICE_DEFAULT)

    z = base.height.astype(F64)
    bump = amp[:, None] * np.exp(-((z[None, :] - zc[:, None]) / width[:, None]) ** 2)
    t = base.temperature.astype(F64)[None, :] + shift[:, None] + bump
    td = base.dewpoint.astype(F64)[None, :] + shift[:, None] + bump
    in_dry = dry[:, None] & (z[None, :] >= dry_bottom[:, None]) & (z[None, :] <= (dry_bottom + dry_depth)[:, None])
    td = np.where(in_dry, t - dry_depression[:, None], td)
    relh = np.where(in_dry, dry_relh[:, None], base.relh.astype(F64)[None, :])

    for i in range(count):
        description = (f"{name}: shift {shift[i]:+.2f} K, nose {amp[i]:+.2f} K at {zc[i]:.0f} m "
                       f"(width {width[i]:.0f} m)")
        if dry[i]:
            description += f", dry layer {dry_bottom[i]:.0f}-{dry_bottom[i] + dry_depth[i]:.0f} m"
        yield make_case("corpus", description, base.pressure, base.height, t[i], td[i], relh[i],
                        str(dsd_names[i]), rime[i], tice[i])


def corpus_cases(ref, base):
    rng = np.random.default_rng(CORPUS_SEED)
    accepted, skipped = [], 0
    for family in CORPUS_FAMILIES:
        count = family[1]
        kept = 0
        for case in corpus_candidates(base, family, 2 * count, rng):
            case, res = settle(ref, case)
            if res.cloud_top_level == 0:
                skipped += 1
                continue
            accepted.append((case, res))
            kept += 1
            if kept == count:
                break
        if kept < count:
            raise RuntimeError(f"corpus family {family[0]}: only {kept} of {count} cases")
    return accepted, skipped


# ---------------------------------------------------------------- coverage

def _preclassifier_path(case, res):
    tw = res.wetbulb[:res.cloud_top_level + 1].astype(F64)
    top, tice = tw[-1], F64(case.tice)
    if tw.max() <= 273.15:
        if res.precip_type == "SN":
            return "all subfreezing -> SN"
        return "all subfreezing, top warmer than Tice -> FZRA" if top >= tice else \
            "all subfreezing, none below Tice in 3 km -> FZRA"
    return "warm top over a cold surface -> FZRA" if res.precip_type == "FZRA" else "all warm -> RA"


def _cloud_top_rule(case, res):
    """Which rule of ptop_slice set the cloud top (labels only; the top itself comes from ptop_slice)."""
    t, td, rh = (a[::-1].astype(F64) for a in (case.temperature, case.dewpoint, 100.0 * case.relh.astype(F64)))
    moist = np.flatnonzero((t - td <= 6) & (rh > 60))
    if res.cloud_top_level == MISSING:
        return "no cloud"
    k = len(t) - 1 - res.cloud_top_level
    if len(moist) == 0:
        return "RH >= 80 % fallback"
    return "dry-layer restart" if k != moist[0] else "first moist level"


def _pl_remelts(res):
    if res.stage != "core":
        return False
    for k, branch in enumerate(res.branch_by_k):
        if branch == "F" and k > 0 and np.any((res.psd_ptype[k - 1] == 4.5) & (res.water_fraction[k - 1] != 1)):
            return True
    return False


def coverage_rows():
    """(item, predicate) pairs that the named group must all satisfy at least once."""
    def core(case, res):
        return res.stage == "core"

    rows = [(f"category {c}", (lambda c: lambda case, res: res.precip_type == c)(c)) for c in PRECIP_TYPE]
    rows += [(f"core category {c}", (lambda c: lambda case, res: core(case, res) and res.precip_type == c)(c))
             for c in PRECIP_TYPE]
    for path in ("all subfreezing -> SN", "all subfreezing, top warmer than Tice -> FZRA",
                 "all subfreezing, none below Tice in 3 km -> FZRA",
                 "warm top over a cold surface -> FZRA", "all warm -> RA"):
        rows.append((f"pre-classifier: {path}", (lambda p: lambda case, res: res.stage == "preclassifier"
                                                 and _preclassifier_path(case, res) == p)(path)))
    rows += [(f"branch {b}", (lambda b: lambda case, res: b in res.branches)(b)) for b in "ABCDEFG"]
    for subset in ("ABF", "ABFG", "CDF"):
        rows.append((f"branches a subset of {{{','.join(subset)}}}",
                     (lambda s: lambda case, res: core(case, res) and set(res.branches) <= set(s))(subset)))
    rows += [
        ("Tice switch", lambda case, res: res.tnuc_switched),
        ("shared refreeze level moved", lambda case, res: res.refreeze_level_moved),
        ("no cloud", lambda case, res: res.stage == "no_cloud"),
        ("cloud top at the highest level", lambda case, res: res.cloud_top_level == len(case.height) - 1),
        ("dry-layer restart", lambda case, res: _cloud_top_rule(case, res) == "dry-layer restart"),
        ("RH >= 80 % fallback", lambda case, res: _cloud_top_rule(case, res) == "RH >= 80 % fallback"),
    ]
    for subset in ("ABF", "ABFG"):
        def within(case, res, s=subset):
            return core(case, res) and set(res.branches) <= set(s)
        label = f"{{{','.join(subset)}}}"
        rows += [(f"{label} with DSD {d}", (lambda d, w=within: lambda case, res: w(case, res)
                                             and case.dsd_name == d)(d)) for d in DSDS]
        rows += [
            (f"{label} with rime factor 5", lambda case, res, w=within: w(case, res) and case.rime_factor == 5),
            (f"{label} with Tice 263.15 K", lambda case, res, w=within: w(case, res)
             and case.tice == F32(TICE_ALT)),
        ]
    rows += [
        ("Nc = 1, warm surface (core RA/SN/RASN)", lambda case, res: core(case, res) and res.crossings == 1
         and res.wetbulb[0] > F32(273.15)),
        ("refrozen PL melts again (A, F, G, F)", lambda case, res: _pl_remelts(res)),
        ("supercooled liquid top over a warm surface", lambda case, res: core(case, res) and "C" in res.branches
         and res.wetbulb[res.cloud_top_level] <= F32(273.15) and res.wetbulb[0] > F32(273.15)),
    ]
    return rows


def print_named_coverage(named):
    print("\nNamed coverage (item, count, case_ids)")
    missing = []
    for item, predicate in coverage_rows():
        ids = [cid for cid, case, res in named if predicate(case, res)]
        print(f"  {item:<64s} {len(ids):3d}  {ids}")
        if not ids:
            missing.append(item)
    if missing:
        raise RuntimeError(f"named cases miss: {missing}")


def print_corpus_coverage(corpus):
    print(f"\nCorpus: {len(corpus)} cases")
    stages = pd.Series([res.stage for _, _, res in corpus]).value_counts()
    print("  stages: " + ", ".join(f"{k} {v}" for k, v in stages.items()))
    for stage in ("preclassifier", "core"):
        counts = {c: sum(1 for _, _, r in corpus if r.stage == stage and r.precip_type == c) for c in PRECIP_TYPE}
        print(f"  {stage:<13s} " + "  ".join(f"{c} {n}" for c, n in counts.items()))
    core = [r for _, _, r in corpus if r.stage == "core"]
    letters = {b: sum(b in r.branches for r in core) for b in "ABCDEFG"}
    print("  branch letters: " + "  ".join(f"{b} {n}" for b, n in letters.items()))
    sets = pd.Series([r.branches for r in core]).value_counts()
    print("  branch sets: " + ", ".join(f"{k} {v}" for k, v in sets.items()))
    print(f"  Tice switch {sum(r.tnuc_switched for r in core)}, "
          f"shared refreeze level moved {sum(r.refreeze_level_moved for r in core)}")
    by_dsd = pd.Series([c.dsd_name for _, c, _ in corpus]).value_counts()
    print("  DSDs: " + ", ".join(f"{k} {v}" for k, v in by_dsd.items())
          + f"; rime 5: {sum(c.rime_factor == 5 for _, c, _ in corpus)}"
          + f"; Tice 263.15 K: {sum(c.tice == F32(TICE_ALT) for _, c, _ in corpus)}")
    lacking = [c for c in PRECIP_TYPE if not any(r.precip_type == c for r in core)]
    lacking += [f"branch {b}" for b, n in letters.items() if n == 0]
    if lacking:
        raise RuntimeError(f"corpus core runs miss: {lacking}")


def print_flags(corpus):
    n = len(corpus)
    near_t = [(cid, r) for cid, _, r in corpus if r.near_threshold]
    near_d = [(cid, r) for cid, _, r in corpus if r.near_discontinuity]
    print(f"\nCorpus flags: near_threshold {len(near_t)} ({100 * len(near_t) / n:.2f} %), "
          f"near_discontinuity {len(near_d)} ({100 * len(near_d) / n:.2f} %)")
    for cid, r in near_d:
        print(f"  near_discontinuity case {cid}: switch {r.disc_switch}, bin {r.disc[0]}, "
              f"level {r.disc[1]}, scope {r.disc[2]}")
    if len(near_d) > CORPUS_CAP * n:
        raise RuntimeError("near_discontinuity exceeds 1 % of the corpus")


def check_deld(ref, rows):
    """Reruns every core case with deld = 1 and requires the same results."""
    worst_lf = worst_fw = 0.0
    n = 0
    for cid, case, res in rows:
        if res.stage != "core":
            continue
        other = run_reference(ref, case, deld=DELD_CHECK)
        same = (other.precip_type, other.supercooled_liquid_height, other.crossings, other.branch_by_k,
                other.near_threshold, other.near_discontinuity, other.disc, other.tnuc_switched,
                other.refreeze_level_moved) == \
               (res.precip_type, res.supercooled_liquid_height, res.crossings, res.branch_by_k,
                res.near_threshold, res.near_discontinuity, res.disc, res.tnuc_switched,
                res.refreeze_level_moved)
        worst_lf = max(worst_lf, abs(other.liquid_fraction - res.liquid_fraction))
        worst_fw = max(worst_fw, float(np.abs(other.water_fraction - res.water_fraction).max()))
        if not same or worst_lf > 1e-12 or worst_fw != 0.0:
            raise RuntimeError(f"case {cid} depends on deld")
        n += 1
    print(f"\ndeld check: {n} core cases rerun with deld = {DELD_CHECK}: identical categories, heights, "
          f"branches, and flags; max |d liquid_fraction| {worst_lf:.3g}, max |d profile| {worst_fw:.3g}")


def print_sample(cid, case, res):
    top = res.cloud_top_level
    print(f"\nSample case {cid}: cloud top level {top}, {case.height[top]:.3f} m AGL, "
          f"{case.pressure[top] / 100:.0f} hPa; {res.precip_type}; liquid fraction {res.liquid_fraction!r}; "
          f"supercooled-liquid height {res.supercooled_liquid_height!r} m AGL; branches {res.branches}; "
          f"crossings {res.crossings}")
    if res.near_threshold:
        print(f"  flagged near_threshold: liquid fraction {res.liquid_fraction}")
    if res.near_discontinuity:
        print(f"  flagged near_discontinuity: switch {res.disc_switch}, bin {res.disc[0]}, "
              f"level {res.disc[1]}, scope {res.disc[2]}")
    expected = (F32(9495.977), F32(27500.0), "PL", 0.0, F32(1130.5469))
    got = (case.height[top], case.pressure[top], res.precip_type, res.liquid_fraction,
           F32(res.supercooled_liquid_height))
    if got != expected:
        raise RuntimeError(f"sample case gives {got}, not {expected}")


# ------------------------------------------------------------------ output

def write_parquet(rows, out_dir):
    levels = {k: [] for k in ("case_id", "level", "pressure", "height", "temperature", "dewpoint",
                              "relh", "wetbulb")}
    cases = {k: [] for k in ("case_id", "group", "description", "dsd_name", "rime_factor",
                             "ice_nucleation_temperature", "cloud_top_height", "cloud_top_level",
                             "precip_type", "liquid_fraction", "supercooled_liquid_height", "stage",
                             "crossings", "branches", "tnuc_switched", "refreeze_level_moved",
                             "near_threshold", "near_discontinuity", "disc_bin", "disc_level",
                             "disc_scope")}
    profiles = {k: [] for k in ("case_id", "level", "bin", "liquid_fraction")}
    for cid, case, res in rows:
        n = len(case.height)
        nbins = len(DSDS[case.dsd_name].diameter)
        levels["case_id"].append(np.full(n, cid, np.int32))
        levels["level"].append(np.arange(n, dtype=np.int32))
        for key, arr in (("pressure", case.pressure), ("height", case.height),
                         ("temperature", case.temperature), ("dewpoint", case.dewpoint),
                         ("relh", case.relh), ("wetbulb", res.wetbulb)):
            levels[key].append(arr)
        top = res.cloud_top_level
        cloud = top != MISSING
        values = (cid, case.group, case.description, case.dsd_name, case.rime_factor, case.tice,
                  case.height[top] if cloud else F32(MISSING), top,
                  PRECIP_TYPE[res.precip_type] if cloud else MISSING, res.liquid_fraction,
                  res.supercooled_liquid_height, res.stage, res.crossings, res.branches,
                  res.tnuc_switched, res.refreeze_level_moved, res.near_threshold,
                  res.near_discontinuity, *res.disc)
        for key, value in zip(cases, values):
            cases[key].append(value)
        lf = np.full((n, nbins), MISSING, F32)
        if res.stage == "core":
            lf[:top + 1] = res.water_fraction[::-1].astype(F32)
        profiles["case_id"].append(np.full(n * nbins, cid, np.int32))
        profiles["level"].append(np.repeat(np.arange(n, dtype=np.int32), nbins))
        profiles["bin"].append(np.tile(np.arange(nbins, dtype=np.int32), n))
        profiles["liquid_fraction"].append(lf.ravel())

    i32, f32, f64, s, b = pa.int32(), pa.float32(), pa.float64(), pa.string(), pa.bool_()
    level_types = dict(case_id=i32, level=i32, pressure=f32, height=f32, temperature=f32,
                       dewpoint=f32, relh=f32, wetbulb=f32)
    case_types = dict(case_id=i32, group=s, description=s, dsd_name=s, rime_factor=f32,
                      ice_nucleation_temperature=f32, cloud_top_height=f32, cloud_top_level=i32,
                      precip_type=i32, liquid_fraction=f64, supercooled_liquid_height=f32, stage=s,
                      crossings=i32, branches=s, tnuc_switched=b, refreeze_level_moved=b,
                      near_threshold=b, near_discontinuity=b, disc_bin=i32, disc_level=i32,
                      disc_scope=s)
    profile_types = dict(case_id=i32, level=i32, bin=i32, liquid_fraction=f32)
    dsd_rows = {"dsd_name": [], "bin": [], "diameter": [], "concentration": []}
    for name, dsd in DSDS.items():
        dsd_rows["dsd_name"] += [name] * len(dsd.diameter)
        dsd_rows["bin"] += list(range(len(dsd.diameter)))
        dsd_rows["diameter"] += list(dsd.diameter)
        dsd_rows["concentration"] += list(dsd.concentration)
    dsd_types = dict(dsd_name=s, bin=i32, diameter=f32, concentration=f32)

    def table(columns, types, concat):
        arrays = [pa.array(np.concatenate(columns[k]) if concat else columns[k], type=t)
                  for k, t in types.items()]
        return pa.Table.from_arrays(arrays, names=list(types))

    os.makedirs(out_dir, exist_ok=True)
    sizes = {}
    for name, tbl in (("levels", table(levels, level_types, True)),
                      ("cases", table(cases, case_types, False)),
                      ("profiles", table(profiles, profile_types, True)),
                      ("dsds", table(dsd_rows, dsd_types, False))):
        path = os.path.join(out_dir, f"{name}.parquet")
        pq.write_table(tbl, path, compression="zstd")
        sizes[name] = (tbl.num_rows, os.path.getsize(path))
    return sizes


def check_dsds():
    default = DSDS["python_default"]
    if not (np.array_equal(default.diameter, np.array([0.05, 0.75, 1.45, 2.15], F32)) and
            np.array_equal(default.concentration, np.array([55.1843, 146.647, 11.6891, 3.60886], F32))):
        raise RuntimeError("python_default does not round to the float32 default DSD literals")
    if len(DSDS["python_deld_0.1"].diameter) != 19:
        raise RuntimeError("python_deld_0.1 does not have 19 bins")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference-dir", default="/workspace/papers/sbc_python")
    parser.add_argument("--output-dir", default=os.path.dirname(os.path.abspath(__file__)))
    args = parser.parse_args()
    start = time.perf_counter()
    warnings.simplefilter("error", RuntimeWarning)
    check_dsds()
    ref = load_reference(args.reference_dir)

    named = []
    for i, case in enumerate(named_cases()):
        settled, res = clear_named(ref, case)
        if settled.description != case.description:
            print(f"named case {i} adjusted: {settled.description}")
        if res.cloud_top_level == 0:
            raise RuntimeError(f"named case {i} has a one-level column")
        named.append((i, settled, res))
    print_named_coverage(named)

    sample, sample_res = settle(ref, sample_case(ref))
    if sample.description != sample_case(ref).description:
        raise RuntimeError("the sample case has exact-threshold inputs")
    print_sample(1000, sample, sample_res)

    accepted, skipped = corpus_cases(ref, sample)
    corpus = [(10000 + i, case, res) for i, (case, res) in enumerate(accepted)]
    nudged = sum("[nudged" in case.description for _, case, _ in corpus)
    print(f"\nCorpus candidates skipped for a one-level column: {skipped}; nudged: {nudged}")
    print_corpus_coverage(corpus)
    print_flags(corpus)

    rows = named + [(1000, sample, sample_res)] + corpus
    for cid, case, res in rows:
        if threshold_violations(case, res):
            raise RuntimeError(f"case {cid} has exact-threshold inputs")
    check_deld(ref, rows)
    sizes = write_parquet(rows, args.output_dir)
    print(f"\nWrote {args.output_dir}:")
    for name, (nrows, size) in sizes.items():
        print(f"  {name}.parquet  {nrows} rows  {size / 1024:.1f} KiB")
    print(f"  total {sum(s for _, s in sizes.values()) / 1024 / 1024:.2f} MiB; "
          f"runtime {time.perf_counter() - start:.1f} s")


if __name__ == "__main__":
    main()
