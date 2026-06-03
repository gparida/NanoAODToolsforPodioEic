#!/usr/bin/env python3
"""
TMVA BDT training: ALP signal vs QED Compton background.

Reads flat ntuples produced by produce_ntuples.py, trains a Boosted Decision
Tree using ROOT TMVA, and writes the trained weights + MVA score tree.

Usage examples
──────────────
  # Train on a single ALP mass hypothesis
  python3 bdt_train.py --signal ma_1.0

  # Combine all signal masses into one inclusive training sample
  python3 bdt_train.py --signal all

  # Limit events (useful for debugging)
  python3 bdt_train.py --signal ma_1.0 --nevts 5000

  # Override output file name
  python3 bdt_train.py --signal ma_2.0 --output my_bdt.root

Outputs
───────
  weights/ALP_BDT_<signal_key>.weights.xml   — trained BDT weights (TMVA XML)
  weights/ALP_BDT_<signal_key>.root          — TMVA output with MVA score tree

The TMVA output ROOT file can be inspected with:
  root -l weights/ALP_BDT_ma_1.0.root
  new TMVA::TMVAGui("weights/ALP_BDT_ma_1.0.root")

Notes
─────
  • Event weights are set to 1.0 (unweighted training).  To apply generator-
    level or luminosity weights add a "weight" branch to the ntuple and pass it
    via DataLoader.SetWeightExpression().
  • The DataLoader normalises signal and background to equal area by default
    (NormMode=NumEvents).  Change in config.py DATALOADER_OPTIONS.
  • Hyperparameters are in config.py BDT_OPTIONS.
"""
import os
import sys
import argparse

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from config import (
    SIGNAL_FILES, BKG_FILES, NTUPLE_DIR, WEIGHT_DIR, PLOT_DIR,
    FEATURES, SPECTATORS,
    BDT_METHOD_NAME, BDT_OPTIONS, DATALOADER_OPTIONS,
)

import ROOT                         # available in eic-shell
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _ntuple_path(key):
    return os.path.join(NTUPLE_DIR, f"{key}_ntuple.root")


def _load_tree(path, tree_name="events"):
    """Open a ROOT file and return its TTree.  Raises if file/tree missing."""
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Ntuple not found: {path}\n"
            f"Run produce_ntuples.py first."
        )
    tf = ROOT.TFile.Open(path, "READ")
    if not tf or tf.IsZombie():
        raise IOError(f"Cannot open: {path}")
    tree = tf.Get(tree_name)
    if not tree:
        raise IOError(f"TTree '{tree_name}' not found in {path}")
    tree.SetDirectory(0)   # detach from file so file can be kept open
    return tf, tree         # caller must keep tf alive


def _build_signal_chain(signal_keys):
    """
    Build a TChain from one or more signal ntuple files.

    Returns (list_of_TFile, TChain).  The TFile handles must stay alive
    for the duration of TMVA training.
    """
    chain  = ROOT.TChain("events")
    tfiles = []
    for key in signal_keys:
        path = _ntuple_path(key)
        if not os.path.isfile(path):
            print(f"  [WARN] Signal ntuple not found, skipping: {path}")
            continue
        tf = ROOT.TFile.Open(path, "READ")
        chain.Add(path)
        tfiles.append(tf)
        print(f"  Signal added: {path}  "
              f"({ROOT.TFile.Open(path, 'READ').Get('events').GetEntries()} evts)")
    return tfiles, chain


def _build_bkg_chain():
    """Build a TChain from all background ntuple files."""
    chain  = ROOT.TChain("events")
    tfiles = []
    for key in BKG_FILES:
        path = _ntuple_path(key)
        if not os.path.isfile(path):
            print(f"  [WARN] Background ntuple not found, skipping: {path}")
            continue
        tf = ROOT.TFile.Open(path, "READ")
        chain.Add(path)
        tfiles.append(tf)
        print(f"  Background added: {path}  "
              f"({ROOT.TFile.Open(path, 'READ').Get('events').GetEntries()} evts)")
    return tfiles, chain


# ── Main training routine ─────────────────────────────────────────────────────

def train(signal_keys, output_root, method_name):
    os.makedirs(WEIGHT_DIR, exist_ok=True)
    os.makedirs(PLOT_DIR,   exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  TMVA BDT Training")
    print(f"  Signal keys : {signal_keys}")
    print(f"  Output      : {output_root}")
    print(f"  Weights dir : {WEIGHT_DIR}")
    print(f"{'='*60}\n")

    # ── Open output file and create TMVA Factory ──────────────────────
    out_tfile = ROOT.TFile.Open(output_root, "RECREATE")

    ROOT.TMVA.Tools.Instance()
    factory = ROOT.TMVA.Factory(
        method_name,
        out_tfile,
        "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P:AnalysisType=Classification",
    )

    # ── DataLoader ────────────────────────────────────────────────────
    loader = ROOT.TMVA.DataLoader("dataset")

    for feat in FEATURES:
        loader.AddVariable(feat, "F")

    for spec in SPECTATORS:
        loader.AddSpectator(spec, "F")

    # ── Load signal and background trees ──────────────────────────────
    print("Loading signal trees …")
    sig_files, sig_chain = _build_signal_chain(signal_keys)

    print("Loading background trees …")
    bkg_files, bkg_chain = _build_bkg_chain()

    if sig_chain.GetEntries() == 0:
        sys.exit("ERROR: No signal events found.  Check ntuple paths.")
    if bkg_chain.GetEntries() == 0:
        sys.exit("ERROR: No background events found.  Check ntuple paths.")

    loader.AddSignalTree(sig_chain,     1.0)
    loader.AddBackgroundTree(bkg_chain, 1.0)

    # Sentinel-event cuts: require DIS kinematics to be filled (x ≠ -999)
    # and at least 1 photon candidate present.
    sig_cut = ROOT.TCut("x > -999 && n_gamma >= 1")
    bkg_cut = ROOT.TCut("x > -999 && n_gamma >= 1")

    loader.PrepareTrainingAndTestTree(sig_cut, bkg_cut, DATALOADER_OPTIONS)

    # ── Book BDT method ───────────────────────────────────────────────
    print(f"\nBooking method: {method_name}")
    print(f"  Options: {BDT_OPTIONS}\n")
    factory.BookMethod(
        loader,
        ROOT.TMVA.Types.kBDT,
        method_name,
        BDT_OPTIONS,
    )

    # ── Train → Test → Evaluate ───────────────────────────────────────
    print("Training …")
    factory.TrainAllMethods()

    print("Testing …")
    factory.TestAllMethods()

    print("Evaluating …")
    factory.EvaluateAllMethods()

    out_tfile.Close()

    # Move weights XML to WEIGHT_DIR with a descriptive name
    default_xml = os.path.join("dataset", "weights", f"{method_name}_{method_name}.weights.xml")
    target_xml  = os.path.join(WEIGHT_DIR, f"{method_name}.weights.xml")
    if os.path.isfile(default_xml):
        import shutil
        shutil.copy2(default_xml, target_xml)
        print(f"\nWeights copied to: {target_xml}")

    print(f"\n[bdt_train] Done.  TMVA output: {output_root}")
    print(f"  Inspect with: root -l {output_root}")
    print(f"  Then in ROOT: new TMVA::TMVAGui(\"{output_root}\")")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Train a TMVA BDT for ALP signal vs QED Compton background."
    )
    parser.add_argument(
        "--signal",
        required=True,
        metavar="KEY_OR_all",
        help=(
            "Signal mass key from config.SIGNAL_FILES (e.g. 'ma_1.0'), "
            "or 'all' to combine all signal masses."
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        metavar="PATH",
        help="Output ROOT file path (default: weights/<method_name>.root).",
    )
    parser.add_argument(
        "--method",
        default=None,
        metavar="NAME",
        help=f"TMVA method name (default: {BDT_METHOD_NAME}_<signal_key>).",
    )
    args = parser.parse_args()

    # Resolve signal keys
    if args.signal == "all":
        signal_keys = list(SIGNAL_FILES.keys())
        sig_tag     = "all"
    else:
        if args.signal not in SIGNAL_FILES:
            sys.exit(
                f"ERROR: Signal key '{args.signal}' not found.\n"
                f"Available: {list(SIGNAL_FILES.keys())}"
            )
        signal_keys = [args.signal]
        sig_tag     = args.signal

    method_name = args.method or f"{BDT_METHOD_NAME}_{sig_tag}"
    output_root = args.output or os.path.join(WEIGHT_DIR, f"{method_name}.root")

    train(signal_keys, output_root, method_name)


if __name__ == "__main__":
    main()
