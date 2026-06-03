# ALP BDT Analysis Framework

Boosted Decision Tree (BDT) framework for discriminating Axion-Like Particle
(ALP → γγ) signal against QED Compton background at the EIC, using ROOT TMVA
and the NanoAODTools/PODIO event loop.

---

## Table of Contents

1. [Overview](#overview)
2. [Directory Layout](#directory-layout)
3. [Physics & Feature Design](#physics--feature-design)
4. [Software Requirements](#software-requirements)
5. [Setting Up Elsewhere](#setting-up-elsewhere)
6. [Step-by-Step Workflow](#step-by-step-workflow)
   - [Step 1 — Configure paths](#step-1--configure-paths)
   - [Step 2 — Produce flat ntuples](#step-2--produce-flat-ntuples)
   - [Step 3 — Train the BDT](#step-3--train-the-bdt)
   - [Step 4 — Diagnostic plots](#step-4--diagnostic-plots)
7. [Configuration Reference](#configuration-reference)
8. [Output Files](#output-files)
9. [Extending the Framework](#extending-the-framework)
10. [Troubleshooting](#troubleshooting)

---

## Overview

```
PODIO reco files (ePIC/EIC format)
        │
        ▼
alp_flat_ntuple.py  ←─ PODIOModule
        │  one flat ROOT TTree per sample
        ▼
produce_ntuples.py  ←─ driver: runs ntuple production for all samples
        │  ntuples/*.root
        ▼
bdt_train.py        ←─ TMVA Factory + DataLoader
        │  weights/*.xml  weights/*.root
        ▼
bdt_diagnostics.py  ←─ matplotlib plots
        │  plots/*.png
        ▼
  Results: input distributions, correlations,
           BDT score, ROC curve, variable importance
```

The framework follows the existing `NanoAODToolsforPodioEic` conventions:
every analysis step is a `PODIOModule` subclass driven by `PODIOPostProcessor`,
so it integrates naturally with the rest of the EIC analysis software.

---

## Directory Layout

```
AnalysiPyForALP/
├── NanoAODToolsforPodioEic/        ← existing framework (unchanged)
└── ALPBDT/
    ├── config.py                   ← all configuration (paths, cuts, features, hyperparams)
    ├── alp_flat_ntuple.py          ← PODIOModule → flat ROOT TTree producer
    ├── produce_ntuples.py          ← driver: ntuple production for all samples
    ├── bdt_train.py                ← TMVA BDT training
    ├── bdt_diagnostics.py          ← diagnostic plots (matplotlib)
    ├── ntuples/                    ← produced flat ntuples (auto-created)
    ├── weights/                    ← TMVA XML weight files + ROOT output (auto-created)
    └── plots/                      ← diagnostic figures (auto-created)
```

Input PODIO files live in:
```
Input_PODIOfiles/Poolfiles_ganesh_May26_2026/
├── ma_0.1_hadded.root     ← ALP signal, m_a = 0.1 GeV
├── ma_0.2_hadded.root
├── ma_0.5_hadded.root
├── ma_1.0_hadded.root
├── ma_2.0_hadded.root
├── ma_5.0_hadded.root
├── ma_10.0_hadded.root
├── ma_20.0_hadded.root
└── qed_compton_hadd.root  ← QED Compton background
```

---

## Physics & Feature Design

### Signal process

ALP (axion-like particle) produced in electron–ion collisions, decaying to two
photons: `e + A → e + a → e + γ + γ`.  Key signatures:
- Two photon clusters in the calorimeter
- Diphoton invariant mass peaks at `m_a`
- Scattered electron carrying the DIS kinematics

### Background process

QED Compton scattering: `e + γ → e + γ` (Bethe-Heitler / virtual photon).
Has one photon and an electron, but no diphoton resonance.

### Feature list

| Feature | Description |
|---------|-------------|
| `x`, `Q2`, `y` | Bjorken x, momentum transfer Q², inelasticity (DIS kinematics, electron method) |
| `e_pt`, `e_eta`, `e_phi`, `e_E`, `e_pz` | Scattered electron kinematics |
| `n_gamma` | Number of photon candidates (neutral reco particles with E > threshold) |
| `g1_E`, `g1_pt`, `g1_eta`, `g1_phi` | Leading photon kinematics |
| `g2_E`, `g2_pt`, `g2_eta`, `g2_phi` | Sub-leading photon kinematics |
| `m_gg` | Di-photon invariant mass ← **primary ALP discriminator** |
| `pt_gg`, `eta_gg`, `dR_gg`, `dphi_gg`, `deta_gg` | Di-photon system kinematics |
| `n_charged`, `n_neutral` | Charged / neutral particle multiplicities |
| `HT`, `HT_gamma_frac` | Scalar HT; fraction from photons |
| `MET_pt`, `MET_phi` | Missing transverse momentum |
| `sum_E` | Total reconstructed energy |
| `dR_eg1`, `dR_eg2` | ΔR between electron and photons |
| `m_egg` | Invariant mass of (e + γ₁ + γ₂) system |

### Object reconstruction

- **Scattered electron**: highest-pT particle in `ReconstructedParticles`
  with |charge| > 0.5, pT > `PT_MIN_ELECTRON`, |η| < `ETA_MAX_ELECTRON`.
- **Photon candidates**: particles in `ReconstructedParticles` with
  |charge| < 0.5, E > `E_MIN_PHOTON`, |η| < `ETA_MAX_PHOTON`, sorted by E.
- **DIS kinematics**: taken directly from `InclusiveKinematicsElectron`.

Events without DIS kinematics or fewer than 1 photon candidate are excluded
from training via a TMVA preselection cut (`x > -999 && n_gamma >= 1`).

---

## Software Requirements

| Package | Version | Source |
|---------|---------|--------|
| Python  | ≥ 3.9   | eic-shell |
| ROOT    | ≥ 6.24 (with PyROOT + TMVA) | eic-shell |
| matplotlib | ≥ 3.5 | `pip install matplotlib` or eic-shell |
| numpy   | ≥ 1.21  | `pip install numpy` or eic-shell |
| NanoAODToolsforPodioEic | — | `../NanoAODToolsforPodioEic/` |

**Everything must be run inside eic-shell** where ROOT and the ePIC software
stack are available.

---

## Setting Up Elsewhere

### 1. Clone the repository structure

```bash
# Your top-level analysis directory
mkdir ~/MyEICAnalysis && cd ~/MyEICAnalysis

# Copy or clone NanoAODToolsforPodioEic
git clone <your-repo> NanoAODToolsforPodioEic

# Copy ALPBDT/
cp -r /path/to/ALPBDT ./ALPBDT
```

### 2. Edit `config.py`

The only file you need to change for a different machine / data location:

```python
# In ALPBDT/config.py — update _PODIO_DIR to your data location
_PODIO_DIR = "/your/path/to/PODIO_files"

# Update NANO_ROOT if NanoAODToolsforPodioEic is in a different location
# (it is computed automatically relative to _HERE by default)
```

Everything else (ntuple output, weights, plots) is written under `ALPBDT/`
itself and requires no path changes.

### 3. Verify NanoAODTools can be imported

```bash
cd ~/MyEICAnalysis/ALPBDT
python3 -c "from alp_flat_ntuple import ALPFlatNtupleProducer; print('OK')"
```

---

## Step-by-Step Workflow

### Step 1 — Configure paths

Edit `config.py` to point `_PODIO_DIR` at your PODIO file directory.
All other paths are derived automatically.

```bash
# (optional) verify all configured signal files exist
python3 - <<'EOF'
import config, os
for k, v in {**config.SIGNAL_FILES, **config.BKG_FILES}.items():
    status = "OK" if os.path.isfile(v) else "MISSING"
    print(f"  [{status}] {k}: {v}")
EOF
```

---

### Step 2 — Produce flat ntuples

This step reads PODIO files and writes a flat ROOT TTree for each sample.
Each tree has one row per event with all features pre-computed.

```bash
# All samples at once (recommended)
python3 produce_ntuples.py --all

# One signal mass only (faster for testing)
python3 produce_ntuples.py --signal ma_1.0

# Quick test: first 2000 events per sample
python3 produce_ntuples.py --all --nevts 2000

# Background only
python3 produce_ntuples.py --bkg
```

Output files (one per sample):
```
ALPBDT/ntuples/ma_0.1_ntuple.root
ALPBDT/ntuples/ma_0.2_ntuple.root
...
ALPBDT/ntuples/qed_compton_ntuple.root
```

Each ROOT file contains a TTree named `events` with branches for every
feature in `config.FEATURES` plus `label` (1=signal, 0=bkg) and
`alp_mass_label` (ALP mass in GeV, 0 for background).

You can inspect a ntuple interactively:
```bash
root -l ntuples/ma_1.0_ntuple.root
# In ROOT prompt:
events->Print()
events->Draw("m_gg")
events->Draw("m_gg>>hgg(100,0,5)", "label==1")
```

---

### Step 3 — Train the BDT

Requires: ntuples from Step 2.

```bash
# Train on one ALP mass hypothesis
python3 bdt_train.py --signal ma_1.0

# Train on all masses combined (inclusive signal)
python3 bdt_train.py --signal all

# Custom output name
python3 bdt_train.py --signal ma_2.0 --output my_ma2_bdt.root
```

Training output:
```
ALPBDT/weights/BDT_ALP_ma_1.0.root       ← TMVA output (MVA scores, trees)
ALPBDT/weights/BDT_ALP_ma_1.0.weights.xml ← trained BDT weights
```

**Interactive TMVA GUI** (requires ROOT with display):
```bash
root -l weights/BDT_ALP_ma_1.0.root
# In ROOT prompt:
new TMVA::TMVAGui("weights/BDT_ALP_ma_1.0.root")
```
The GUI shows the TMVA-standard plots (overtraining check, input variables,
ROC, output distributions).

**Apply weights to new data** (scoring):
```python
import ROOT
reader = ROOT.TMVA.Reader()
# Add variables in the same order as training
reader.AddVariable("x",    x_buf)
reader.AddVariable("m_gg", mgg_buf)
# ...
reader.BookMVA("BDT_ALP_ma_1.0", "weights/BDT_ALP_ma_1.0.weights.xml")
score = reader.EvaluateMVA("BDT_ALP_ma_1.0")
```

---

### Step 4 — Diagnostic plots

Produces matplotlib figures in `ALPBDT/plots/`.

```bash
# Full diagnostics (needs TMVA output from Step 3)
python3 bdt_diagnostics.py --signal ma_1.0

# Only input distributions + correlations (no TMVA file needed)
python3 bdt_diagnostics.py --signal ma_1.0 --correlations-only

# Limit events read (faster)
python3 bdt_diagnostics.py --signal ma_1.0 --nevts 10000
```

Produced plots:

| File | Content |
|------|---------|
| `input_variables_<key>.png` | All feature distributions, signal vs background overlay |
| `correlations_<key>.png` | Pearson ρ correlation matrices for signal and background |
| `bdt_score_<key>.png` | BDT response distributions (signal / background, test sample) |
| `roc_<key>.png` | ROC curve with AUC |
| `variable_importance_<key>.png` | TMVA variable importance ranking |

---

## Configuration Reference

All tunable parameters live in `config.py`.

### File paths
```python
_PODIO_DIR    # directory containing PODIO ROOT files
NTUPLE_DIR    # where flat ntuples are written
WEIGHT_DIR    # where TMVA weights are written
PLOT_DIR      # where diagnostic plots are written
```

### Collection names
```python
COLL_RECON_ALL     = "ReconstructedParticles"
COLL_RECON_CHARGED = "ReconstructedChargedParticles"
COLL_KIN_ELECTRON  = "InclusiveKinematicsElectron"
```
Update these if the ePIC/EIC software version changes the collection naming.

### Object selection thresholds
```python
ETA_MAX_ELECTRON = 4.0   # |η| acceptance for scattered electron
PT_MIN_ELECTRON  = 0.5   # pT threshold (GeV)
E_MIN_PHOTON     = 0.2   # minimum photon energy (GeV)
ETA_MAX_PHOTON   = 4.0   # |η| acceptance for photons
```

### BDT training features
`FEATURES` is a Python list of branch names.  Comment out any entry to
exclude it from BDT training without touching any other script.

### BDT hyperparameters
```python
BDT_OPTIONS = (
    "!H:!V"
    ":NTrees=500"         # number of trees in the ensemble
    ":MaxDepth=3"         # tree depth (controls complexity)
    ":MinNodeSize=2.5%"   # minimum fraction of events per leaf
    ":BoostType=AdaBoost" # boosting algorithm
    ":AdaBoostBeta=0.5"   # learning rate
    ":UseBaggedBoost=True"
    ":BaggedSampleFraction=0.5"
    ":SeparationType=GiniIndex"
    ":nCuts=20"           # grid points for cut optimisation
    ":PruneMethod=NoPruning"
)
```

Useful starting points for tuning:
- Increase `NTrees` (800–1000) and decrease `AdaBoostBeta` (0.1–0.3) for
  lower overtraining risk.
- Decrease `MaxDepth` to 2 if the overtraining check shows large train/test
  separation.
- Use `BoostType=GradBoost` with `Shrinkage=0.1` as an alternative to AdaBoost.

---

## Output Files

### Flat ntuples (`ntuples/*.root`)

One file per sample.  TTree name: `events`.  Key branches:

| Branch | Type | Description |
|--------|------|-------------|
| `label` | int | 1 = ALP signal, 0 = QED Compton background |
| `alp_mass_label` | float | ALP mass hypothesis (GeV); 0 for background |
| `x`, `Q2`, `y` | float | DIS kinematics |
| `m_gg` | float | Di-photon invariant mass (GeV) |
| `n_gamma` | int | Photon candidate multiplicity |
| *(+ all features in `config.FEATURES`)* | float | See feature table above |

Sentinel value for missing / undefined quantities: **-999.0**.

### TMVA weights (`weights/*.xml`)

Standard TMVA XML format.  Can be loaded with `TMVA::Reader` to score new
events.  Compatible with ROOT ≥ 6.

### TMVA output (`weights/*.root`)

Contains TMVA-internal trees for the training / test samples with the
per-event MVA score.  Open with `TMVA::TMVAGui` for the full suite of
TMVA-standard diagnostic plots.

---

## Extending the Framework

### Add a new feature

1. Add the feature computation to `alp_flat_ntuple.py` (`_fill_particles` or
   `_fill_dis_kin`).  Extend `_BRANCHES` at the top of the file.
2. Add the branch name to `config.FEATURES`.
3. Re-run `produce_ntuples.py --all` to regenerate the ntuples.
4. Re-run `bdt_train.py` to retrain.

### Add a new signal sample

1. Add an entry to `config.SIGNAL_FILES`:
   ```python
   "ma_50.0": os.path.join(_PODIO_DIR, "ma_50.0_hadded.root"),
   ```
2. Run `produce_ntuples.py --signal ma_50.0`.
3. Run `bdt_train.py --signal ma_50.0`.

### Train mass-decorrelated BDT

To train a BDT that is insensitive to the diphoton mass (useful for
model-independent searches):
- Add `m_gg` to `config.SPECTATORS` (not `FEATURES`).
- Apply the `Decorrelate` transform in `BDT_OPTIONS`:
  `":Transformations=D"` where `D` stands for decorrelation.

### Add event weights

1. Add a `weight` branch to `alp_flat_ntuple.py`.
2. In `bdt_train.py` add before `PrepareTrainingAndTestTree`:
   ```python
   loader.SetSignalWeightExpression("weight")
   loader.SetBackgroundWeightExpression("weight")
   ```

---

## Troubleshooting

**`ModuleNotFoundError: No module named 'PhysicsTools'`**
- Make sure `NanoAODToolsforPodioEic/` is at the path expected by `config.NANO_ROOT`.
- Check `config.py: NANO_ROOT` is correct.

**`IOError: Cannot open file: ... (zombie)`**
- The PODIO file path in `config.py` is wrong.  Run the path-check snippet in
  Step 1 to find missing files.

**`TTree 'events' not found`**
- The ntuple file exists but is empty / corrupt.  Delete it and rerun
  `produce_ntuples.py`.

**TMVA training: `ERROR: No signal events found`**
- The ntuple file for the signal sample does not exist.  Run
  `produce_ntuples.py --signal <key>` first.

**TMVA training: `Preselection: all events cut`**
- All events have `x == -999`, meaning `InclusiveKinematicsElectron` was
  never filled.  Check the PODIO collection name in `config.COLL_KIN_ELECTRON`
  matches what is in your files:
  ```bash
  python3 ../NanoAODToolsforPodioEic/scripts/run_podio_example.py \
      --input /path/to/file.root --list
  ```

**Plots are blank / all events have sentinel values**
- Check that `config.COLL_RECON_ALL` matches the reconstructed-particle
  collection name in your PODIO files (see `--list` above).

**Overtraining observed in TMVA GUI**
- Reduce `MaxDepth` (try 2) or increase `MinNodeSize` (try 5%).
- Increase `NTrees` and reduce `AdaBoostBeta` (try 0.2).
- Ensure signal and background samples have similar statistics.

---

*Framework written for the EIC ALP → γγ analysis.  ROOT TMVA reference:
https://root.cern/manual/tmva/*
