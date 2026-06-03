#!/usr/bin/env python3
"""
Produce flat ROOT ntuples from PODIO signal and background files.

Runs ALPFlatNtupleProducer on every sample configured in config.py and
writes one flat ROOT file per sample to ALPBDT/ntuples/.

Usage examples
──────────────
  # Produce ntuples for ALL samples (signal + background)
  python3 produce_ntuples.py --all

  # One specific signal mass
  python3 produce_ntuples.py --signal ma_1.0

  # Background only
  python3 produce_ntuples.py --bkg

  # Limit events per sample (useful for testing)
  python3 produce_ntuples.py --all --nevts 2000

Output
──────
  ntuples/ma_0.1_ntuple.root
  ntuples/ma_0.2_ntuple.root
  ...
  ntuples/qed_compton_ntuple.root

Each file contains one TTree named "events" with the branches defined
in alp_flat_ntuple.py.  The "label" branch is 1 for signal, 0 for bkg.
"""
import os
import sys
import argparse

# ── Ensure this directory is on sys.path ──────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

# alp_flat_ntuple bootstraps NanoAODTools; config just needs os/sys.
from config import (
    SIGNAL_FILES, BKG_FILES, NTUPLE_DIR,
    COLL_RECON_ALL, COLL_KIN_ELECTRON,
    ETA_MAX_ELECTRON, PT_MIN_ELECTRON,
    E_MIN_PHOTON, ETA_MAX_PHOTON,
    NANO_ROOT,
)
from alp_flat_ntuple import ALPFlatNtupleProducer


def _add_nanoaod_path():
    """Register PODIOPostProcessor on the import path."""
    import types as _types
    _pyroot  = os.path.join(NANO_ROOT, "python")
    _physdir = os.path.join(NANO_ROOT, "build", "lib", "python", "PhysicsTools")

    def _ns(name, path):
        if name not in sys.modules:
            m = _types.ModuleType(name)
            m.__path__ = [path]
            m.__package__ = name
            sys.modules[name] = m

    _ns("PhysicsTools",              _physdir)
    _ns("PhysicsTools.NanoAODTools", _pyroot)
    if _pyroot not in sys.path:
        sys.path.insert(0, _pyroot)


def run_sample(name, podio_file, label, alp_mass, nevts):
    """Run the ntuple producer for a single PODIO file."""
    if not os.path.isfile(podio_file):
        print(f"  [WARN] File not found, skipping: {podio_file}")
        return None

    _add_nanoaod_path()
    from PhysicsTools.NanoAODTools.postprocessing.framework.podio_postprocessor import PODIOPostProcessor

    out_file = os.path.join(NTUPLE_DIR, f"{name}_ntuple.root")
    print(f"\n{'='*60}")
    print(f"  Sample  : {name}")
    print(f"  Input   : {podio_file}")
    print(f"  Output  : {out_file}")
    print(f"  Label   : {label}  ALP mass : {alp_mass} GeV")
    print(f"  Max evts: {nevts if nevts else 'all'}")
    print(f"{'='*60}")

    producer = ALPFlatNtupleProducer(
        output_file   = out_file,
        label         = label,
        alp_mass      = alp_mass,
        coll_recon    = COLL_RECON_ALL,
        coll_kin      = COLL_KIN_ELECTRON,
        eta_max_e     = ETA_MAX_ELECTRON,
        pt_min_e      = PT_MIN_ELECTRON,
        e_min_gamma   = E_MIN_PHOTON,
        eta_max_gamma = ETA_MAX_PHOTON,
    )

    PODIOPostProcessor(
        inputFiles    = [podio_file],
        modules       = [producer],
        maxEntries    = nevts,
        progressEvery = 500,
    ).run()

    return out_file


def main():
    parser = argparse.ArgumentParser(
        description="Produce flat ROOT ntuples from PODIO signal/background files."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--all",
        action="store_true",
        help="Process all signal masses and the background sample.",
    )
    group.add_argument(
        "--signal",
        metavar="MASS_KEY",
        help="Process one signal sample by key, e.g. ma_1.0",
    )
    group.add_argument(
        "--bkg",
        action="store_true",
        help="Process only the background sample(s).",
    )
    parser.add_argument(
        "--nevts",
        type=int,
        default=None,
        metavar="N",
        help="Maximum events to process per file (default: all).",
    )
    args = parser.parse_args()

    os.makedirs(NTUPLE_DIR, exist_ok=True)

    # Build the work list: list of (name, podio_file, label, alp_mass)
    tasks = []

    if args.all or args.signal:
        if args.signal:
            if args.signal not in SIGNAL_FILES:
                sys.exit(
                    f"ERROR: Signal key '{args.signal}' not found.\n"
                    f"Available: {list(SIGNAL_FILES.keys())}"
                )
            keys = [args.signal]
        else:
            keys = list(SIGNAL_FILES.keys())

        for key in keys:
            # Extract the numeric mass from the key name (e.g. "ma_1.0" → 1.0)
            try:
                mass = float(key.split("_", 1)[1])
            except (IndexError, ValueError):
                mass = 0.0
            tasks.append((key, SIGNAL_FILES[key], 1, mass))

    if args.all or args.bkg:
        for key, path in BKG_FILES.items():
            tasks.append((key, path, 0, 0.0))

    produced = []
    for name, path, label, mass in tasks:
        out = run_sample(name, path, label, mass, args.nevts)
        if out:
            produced.append(out)

    print(f"\n[produce_ntuples] Done. {len(produced)} file(s) written to {NTUPLE_DIR}/")
    for f in produced:
        print(f"  {f}")


if __name__ == "__main__":
    main()
