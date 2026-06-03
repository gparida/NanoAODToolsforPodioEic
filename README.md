# NanoAODTools for EIC/ePIC PODIO Files

A standalone Python analysis framework for reading EIC/ePIC reconstruction files in PODIO format,
built on top of the NanoAODTools infrastructure.
No CMSSW installation is required — only ROOT with PyROOT bindings.

---

## What this is

This repo adapts the NanoAODTools event-loop and module system to work directly on PODIO ROOT files
produced by the ePIC simulation/reconstruction chain (e.g. `recon_*.root`).
The API mirrors the NanoAODTools style (`Collection`, `Module`, `PostProcessor`) so that analysis code
is easy to read and port.

Key features:
- Read any PODIO `events` tree without loading EDM4hep/PODIO C++ dictionaries
- Per-event collection access with attribute-style member lookup (`p.energy`, `p.momentum_x`)
- Built-in helpers: `pt()`, `eta()`, `phi()`, `p4()`
- Rec↔MC truth matching via `RecoMCAssociation`
- MC parent/daughter navigation via `MCParticleNavigator`

---

## Requirements

| Dependency | Where to get it |
|---|---|
| Python 3.6+ | Included in eic-shell |
| ROOT + PyROOT | Included in eic-shell |
| matplotlib *(optional)* | `pip install matplotlib` — only needed for plot output |

**eic-shell** (the standard EIC software container) provides everything mandatory.
No `pip install` is needed for the core framework.

---

## Setup on a Remote EIC Cluster

These instructions apply to BNL SDCC/RCF, JLab, NERSC, or any cluster that provides
the EIC software container via Singularity/Apptainer.

### Step 1 — Log in and clone the repository

```bash
ssh <your_username>@<cluster_hostname>

# Clone into your home or work directory — the repo is lightweight (~few MB, no large files)
git clone https://github.com/gparida/NanoAODToolsforPodioEic.git
cd NanoAODToolsforPodioEic
```

Recommended locations:
- **Scripts/code** → `$HOME/` or `$HOME/analysis/` (backed up, small quota is fine)
- **Input recon files** → `/gpfs/...`, `/lustre/...`, or `/work/...` (large storage partition)
- **Output files** → same large-storage partition as inputs, or a dedicated `output/` directory

### Step 2 — Enter the EIC software environment

The EIC software container (eic-shell) provides ROOT, PyROOT, and all physics libraries.

**On BNL SDCC/RCF:**
```bash
# Option A: interactive shell
singularity shell /cvmfs/eic.opensciencegrid.org/singularity/epic-simulation-environment_latest.sif

# Option B: run a single command inside the container
singularity exec /cvmfs/eic.opensciencegrid.org/singularity/epic-simulation-environment_latest.sif \
    python scripts/run_podio_example.py --input /path/to/recon.root
```

**On JLab / NERSC / generic cluster with Apptainer:**
```bash
apptainer shell /path/to/eic-shell.sif
# or use the eic-shell alias if your site has set it up:
eic-shell
```

Once inside the container, verify ROOT is available:
```bash
python -c "import ROOT; print(ROOT.__version__)"
```

### Step 3 — Run the example analysis

```bash
# Inside eic-shell, from the repo root directory:
python scripts/run_podio_example.py --input /path/to/recon.root

# Limit to first 100 events (useful for a quick test):
python scripts/run_podio_example.py --input /path/to/recon.root --nevts 100

# Skip the first 500 events:
python scripts/run_podio_example.py --input /path/to/recon.root --first 500

# List all collections available in the file (no analysis, just inspection):
python scripts/run_podio_example.py --input /path/to/recon.root --list
```

Expected output (truncated):
```
Input file : /path/to/recon.root
[PODIOPostProcessor] Opening: /path/to/recon.root
  990 events selected (firstEntry=0, maxEntries=None)
  Collections available: 42
  Processed     100/990 events  (10.1%)  rate=12.3 kHz  accepted=100
  ...
Done. Processed 990 events in 0.8s (1.24 kHz). Accepted 990/990.

[EICExampleAnalysis] Results after 990 events:
  ReconstructedParticles per event : mean=18.42  min=0  max=61
  Recon |p|  : mean=1.847 GeV
  ...
```

### Step 4 — Run as a batch job

Create a submission script (example for Slurm):

```bash
#!/bin/bash
#SBATCH --job-name=eic_analysis
#SBATCH --output=logs/eic_%j.out
#SBATCH --time=02:00:00
#SBATCH --mem=4G

singularity exec /cvmfs/eic.opensciencegrid.org/singularity/epic-simulation-environment_latest.sif \
    python /path/to/NanoAODToolsforPodioEic/scripts/run_podio_example.py \
        --input /gpfs/mydata/recon_170.root \
        --nevts 10000
```

Submit with:
```bash
sbatch my_job.sh
```

---

## Writing Your Own Analysis

### Step 1 — Create a new module file

Create a file, e.g. `python/postprocessing/examples/my_analysis.py`:

```python
import math
from PhysicsTools.NanoAODTools.postprocessing.framework.podio_eventloop import PODIOModule
from PhysicsTools.NanoAODTools.postprocessing.framework.podio_datamodel import Collection


class MyAnalysis(PODIOModule):

    def beginJob(self):
        self.electron_pt = []

    def analyze(self, event):
        charged = Collection(event, "ReconstructedChargedParticles")
        for p in charged:
            if abs(p.energy) > 0.1:                     # basic energy cut
                self.electron_pt.append(p.pt())
        return True                                      # keep all events

    def endJob(self):
        if self.electron_pt:
            print("Mean pT of charged particles: {:.3f} GeV".format(
                sum(self.electron_pt) / len(self.electron_pt)
            ))
```

### Step 2 — Create a runner script

Create `scripts/run_my_analysis.py`:

```python
#!/usr/bin/env python3
import os, sys, types, argparse

# --- path bootstrap (copy this block into every new runner script) ---
_here         = os.path.dirname(os.path.abspath(__file__))
_repo_root    = os.path.dirname(_here)
_python_src   = os.path.join(_repo_root, "python")

def _ns(name, path):
    if name not in sys.modules:
        m = types.ModuleType(name); m.__path__ = [path]; sys.modules[name] = m
_ns("PhysicsTools",              os.path.join(_repo_root, "build", "lib", "python", "PhysicsTools"))
_ns("PhysicsTools.NanoAODTools", _python_src)
if _python_src not in sys.path:
    sys.path.insert(0, _python_src)
# ---

from PhysicsTools.NanoAODTools.postprocessing.framework.podio_postprocessor import PODIOPostProcessor
from PhysicsTools.NanoAODTools.postprocessing.examples.my_analysis import MyAnalysis

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--nevts", type=int, default=None)
args = parser.parse_args()

PODIOPostProcessor(
    inputFiles=[args.input],
    modules=[MyAnalysis()],
    maxEntries=args.nevts,
).run()
```

Run it:
```bash
python scripts/run_my_analysis.py --input /path/to/recon.root
```

### Available collections (typical ePIC recon file)

| Collection name | Contents |
|---|---|
| `ReconstructedParticles` | All reconstructed particles (charged + neutral) |
| `ReconstructedChargedParticles` | Charged-only reconstructed tracks |
| `MCParticles` | Generator-level truth particles |
| `InclusiveKinematicsElectron` | DIS kinematics via electron method (x, Q2, y) |
| `ReconstructedChargedParticleAssociations` | Rec↔MC truth links |

Common members per particle:
- `energy`, `momentum.x`, `momentum.y`, `momentum.z`, `mass`, `charge`
- `PDG` (MCParticles only), `generatorStatus` (MCParticles only)

Use `--list` to see all collections and their members for a specific file:
```bash
python scripts/run_podio_example.py --input recon.root --list
```

### MC truth matching and parent navigation

```python
from PhysicsTools.NanoAODTools.postprocessing.framework.podio_associations import (
    RecoMCAssociation, MCParticleNavigator
)

def analyze(self, event):
    charged = Collection(event, "ReconstructedChargedParticles")
    mc      = Collection(event, "MCParticles")

    # Rec <-> MC matching
    assoc = RecoMCAssociation(event, "ReconstructedChargedParticleAssociations")
    for i, rec in enumerate(charged):
        mc_idx, weight = assoc.best_mc(i)   # (None, 0.0) if unmatched
        if mc_idx is not None:
            print("PDG =", int(mc[mc_idx]["PDG"]))

    # MC parent/daughter navigation
    nav = MCParticleNavigator(event)
    for i in range(len(mc)):
        if nav.is_stable(i):
            for parent in nav.parents(i):
                print("parent PDG =", int(parent["PDG"]))
```

---

## ALP BDT Analysis

`ALPBDT/` contains a complete TMVA Boosted Decision Tree framework for
discriminating **ALP (axion-like particle) → γγ signal** against
**QED Compton background** at the EIC.

### Quick start

```bash
# 1. Produce flat ROOT ntuples from all PODIO signal + background files
python3 ALPBDT/produce_ntuples.py --all

# 2. Train a BDT for one ALP mass hypothesis
python3 ALPBDT/bdt_train.py --signal ma_1.0

# 3. Diagnostic plots: input distributions, correlations, ROC, variable importance
python3 ALPBDT/bdt_diagnostics.py --signal ma_1.0
```

### Supported signal samples

ALP masses: 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0 GeV.
Background: QED Compton (`qed_compton_hadd.root`).

### BDT input features (32 total)

DIS kinematics (x, Q², y) · scattered-electron kinematics ·
photon multiplicity · leading/sub-leading photon kinematics ·
di-photon system (m_gg, ΔR, Δφ, pT) · event shape (HT, MET, ΣE) ·
electron–photon angular correlations (ΔR, 3-body mass).

### ALPBDT directory layout

```
ALPBDT/
├── config.py              # all paths, cuts, feature list, BDT hyperparameters
├── alp_flat_ntuple.py     # PODIOModule → flat ROOT TTree (one row per event)
├── produce_ntuples.py     # driver: ntuple production for all samples
├── bdt_train.py           # TMVA Factory / DataLoader BDT training
├── bdt_diagnostics.py     # matplotlib: distributions, correlations, ROC, importance
├── README.md              # full setup and usage documentation
├── ntuples/               # produced flat ntuples (git-ignored)
├── weights/               # TMVA XML weight files + ROOT output (git-ignored)
└── plots/                 # diagnostic figures (git-ignored)
```

See **[ALPBDT/README.md](ALPBDT/README.md)** for the full documentation:
configuration reference, step-by-step workflow, feature descriptions,
BDT hyperparameter tuning guide, and troubleshooting.

---

## Repository layout

```
NanoAODToolsforPodioEic/
├── ALPBDT/                          # ALP BDT analysis (see ALPBDT/README.md)
├── python/postprocessing/
│   ├── framework/
│   │   ├── podio_reader.py          # opens ROOT file, iterates events (PyROOT)
│   │   ├── podio_datamodel.py       # PODIOEvent / PODIOCollection / PODIOObject
│   │   ├── podio_eventloop.py       # PODIOModule base class + event loop
│   │   ├── podio_postprocessor.py   # high-level runner
│   │   └── podio_associations.py    # RecoMCAssociation + MCParticleNavigator
│   └── examples/
│       ├── podio_example_analysis.py  # EICExampleAnalysis (used by run_podio_example.py)
│       └── podio_exampleAnalysis.py   # extended example with associations + MC navigation
└── scripts/
    └── run_podio_example.py         # CLI entry point
```

---

## Note on `standalone/env_standalone.sh`

This script is inherited from the original CMS NanoAODTools and is **not needed** for EIC/PODIO
analysis. It requires Python 2.7 and a CMSSW-style build, neither of which is relevant here.
The path bootstrap in each runner script (`scripts/run_podio_example.py`) replaces it entirely.
