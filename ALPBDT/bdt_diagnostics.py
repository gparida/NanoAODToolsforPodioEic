#!/usr/bin/env python3
"""
Post-training diagnostics for the ALP BDT analysis.

Produces matplotlib figures (saved to ALPBDT/plots/):

  1. Input variable distributions — signal vs background overlays (one panel
     per feature, normalised to unit area).

  2. Correlation matrices — Pearson ρ heatmaps for signal and background
     separately.  Help identify redundant features and spot unexpected
     correlations.

  3. BDT score distributions — overlaid signal/background after training,
     read from the TMVA output ROOT file.

  4. ROC curve — background rejection vs signal efficiency, with AUC label.

  5. Variable importance — bar chart of the TMVA-reported variable ranking
     (read from the XML weight file).

Usage examples
──────────────
  # Full diagnostics for one training (after bdt_train.py finishes)
  python3 bdt_diagnostics.py --signal ma_1.0

  # Only the correlation matrices (no TMVA output needed)
  python3 bdt_diagnostics.py --signal ma_1.0 --correlations-only

  # Use a non-default TMVA output file
  python3 bdt_diagnostics.py --signal ma_1.0 --tmva-output my_bdt.root

Notes
─────
  Matplotlib and numpy are used for all plots (no dependency on TMVA's
  built-in GUI, which requires an interactive session).

  ROOT is still needed to read the flat ntuples and the TMVA output file.
"""
import os
import sys
import argparse
import math

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from config import (
    SIGNAL_FILES, BKG_FILES, NTUPLE_DIR, WEIGHT_DIR, PLOT_DIR,
    FEATURES, BDT_METHOD_NAME,
)

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False


# ── TTree → numpy array helpers ───────────────────────────────────────────────

def _tree_to_arrays(path, branches, cut="", max_events=-1):
    """
    Read selected branches from a flat ntuple TTree into a 2-D numpy array.

    Returns
    -------
    data : np.ndarray, shape (n_events, n_branches)
    """
    tf = ROOT.TFile.Open(path, "READ")
    if not tf or tf.IsZombie():
        raise IOError(f"Cannot open: {path}")
    tree = tf.Get("events")
    if not tree:
        raise IOError(f"TTree 'events' not found in {path}")

    # Use ROOT's Draw with a temporary tree to extract arrays
    rows = []
    n = tree.GetEntries()
    if max_events > 0:
        n = min(n, max_events)

    for i, entry in enumerate(tree):
        if i >= n:
            break
        if cut:
            # Simple sentinel cut applied manually
            try:
                x = getattr(entry, "x", -999.0)
                ng = getattr(entry, "n_gamma", 0)
                if x <= -998.0 or ng < 1:
                    continue
            except Exception:
                pass
        row = []
        for b in branches:
            try:
                row.append(float(getattr(entry, b)))
            except Exception:
                row.append(float("nan"))
        rows.append(row)

    tf.Close()
    return np.array(rows, dtype=np.float32) if rows else np.zeros((0, len(branches)), dtype=np.float32)


def _load_sig_bkg(signal_key, max_events=-1):
    """
    Load signal and background arrays from flat ntuples.

    Returns (sig_arr, bkg_arr) both shape (n, n_features).
    """
    sig_path = os.path.join(NTUPLE_DIR, f"{signal_key}_ntuple.root")
    bkg_path = os.path.join(NTUPLE_DIR, "qed_compton_ntuple.root")

    cut = "x > -999 && n_gamma >= 1"

    print(f"  Reading signal  : {sig_path}")
    sig = _tree_to_arrays(sig_path, FEATURES, cut=cut, max_events=max_events)

    print(f"  Reading bkg     : {bkg_path}")
    bkg = _tree_to_arrays(bkg_path, FEATURES, cut=cut, max_events=max_events)

    print(f"  Signal events   : {len(sig)}")
    print(f"  Background events: {len(bkg)}")
    return sig, bkg


# ── Plot 1: Input variable distributions ─────────────────────────────────────

def plot_input_variables(sig, bkg, signal_key, outdir):
    if not _HAS_MPL:
        print("  [SKIP] matplotlib not available")
        return

    n_feat   = len(FEATURES)
    n_cols   = 5
    n_rows   = math.ceil(n_feat / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3.5, n_rows * 3))
    axes = axes.flatten()

    for i, feat in enumerate(FEATURES):
        ax = axes[i]
        s_vals = sig[:, i]
        b_vals = bkg[:, i]

        # Drop sentinel / NaN
        s_vals = s_vals[(s_vals > -998) & np.isfinite(s_vals)]
        b_vals = b_vals[(b_vals > -998) & np.isfinite(b_vals)]

        if len(s_vals) == 0 and len(b_vals) == 0:
            ax.set_title(feat + " (no data)")
            continue

        lo = min(np.percentile(s_vals, 1) if len(s_vals) else 0,
                 np.percentile(b_vals, 1) if len(b_vals) else 0)
        hi = max(np.percentile(s_vals, 99) if len(s_vals) else 1,
                 np.percentile(b_vals, 99) if len(b_vals) else 1)

        bins = np.linspace(lo, hi, 51)
        ax.hist(b_vals, bins=bins, density=True, alpha=0.55,
                color="steelblue", label="QED Compton (bkg)")
        ax.hist(s_vals, bins=bins, density=True, alpha=0.55,
                color="tomato", label=f"ALP {signal_key} (sig)")
        ax.set_xlabel(feat, fontsize=8)
        ax.set_ylabel("Norm. events", fontsize=7)
        ax.tick_params(labelsize=7)
        if i == 0:
            ax.legend(fontsize=7)

    # Hide unused panels
    for j in range(n_feat, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(f"Input variable distributions — {signal_key}", fontsize=11)
    plt.tight_layout()
    out = os.path.join(outdir, f"input_variables_{signal_key}.png")
    plt.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Plot 2: Correlation matrices ──────────────────────────────────────────────

def plot_correlations(sig, bkg, signal_key, outdir):
    if not _HAS_MPL:
        print("  [SKIP] matplotlib not available")
        return

    def _corr_matrix(arr):
        # Replace sentinels with NaN and compute column-wise Pearson ρ
        a = arr.copy().astype(np.float64)
        a[a <= -998] = np.nan
        n = a.shape[1]
        C = np.full((n, n), np.nan)
        for i in range(n):
            for j in range(n):
                mask = np.isfinite(a[:, i]) & np.isfinite(a[:, j])
                if mask.sum() < 2:
                    continue
                C[i, j] = np.corrcoef(a[mask, i], a[mask, j])[0, 1]
        return C

    fig, axes = plt.subplots(1, 2, figsize=(22, 10))

    for ax, arr, title in [
        (axes[0], sig, f"Signal ({signal_key})"),
        (axes[1], bkg, "Background (QED Compton)"),
    ]:
        C = _corr_matrix(arr)
        im = ax.imshow(C, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
        ax.set_xticks(range(len(FEATURES)))
        ax.set_yticks(range(len(FEATURES)))
        ax.set_xticklabels(FEATURES, rotation=90, fontsize=7)
        ax.set_yticklabels(FEATURES, fontsize=7)
        ax.set_title(f"Pearson Correlation — {title}", fontsize=10)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        # Annotate cells
        for i in range(len(FEATURES)):
            for j in range(len(FEATURES)):
                v = C[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                            fontsize=5,
                            color="white" if abs(v) > 0.65 else "black")

    plt.tight_layout()
    out = os.path.join(outdir, f"correlations_{signal_key}.png")
    plt.savefig(out, dpi=120)
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Plot 3 + 4: BDT score and ROC from TMVA output ───────────────────────────

def _read_tmva_scores(tmva_root, method_name):
    """
    Extract signal and background BDT scores from a TMVA output file.

    TMVA stores a TestTree in dataset/TestTree with branch BDT_<name>.
    Returns (sig_scores, bkg_scores) as numpy arrays.
    """
    tf = ROOT.TFile.Open(tmva_root, "READ")
    if not tf or tf.IsZombie():
        raise IOError(f"Cannot open TMVA file: {tmva_root}")

    test_tree = tf.Get("dataset/TestTree")
    if not test_tree:
        # Try alternative location
        test_tree = tf.Get("TestTree")
    if not test_tree:
        tf.Close()
        raise IOError("TestTree not found in TMVA output file.")

    branch_name = method_name   # TMVA stores score under the method name

    sig_scores, bkg_scores = [], []
    for entry in test_tree:
        try:
            score = getattr(entry, branch_name)
            label = int(entry.classID)   # 0=bkg, 1=sig  (TMVA convention)
        except AttributeError:
            continue
        if label == 1:
            sig_scores.append(score)
        else:
            bkg_scores.append(score)

    tf.Close()
    return np.array(sig_scores), np.array(bkg_scores)


def plot_bdt_score(sig_scores, bkg_scores, signal_key, method_name, outdir):
    if not _HAS_MPL:
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    bins = np.linspace(-1, 1, 51)
    ax.hist(bkg_scores, bins=bins, density=True, alpha=0.55,
            color="steelblue", label="QED Compton (bkg)")
    ax.hist(sig_scores, bins=bins, density=True, alpha=0.55,
            color="tomato", label=f"ALP {signal_key} (sig)")
    ax.set_xlabel("BDT response", fontsize=11)
    ax.set_ylabel("Normalised events / bin", fontsize=11)
    ax.set_title(f"BDT score distribution — {signal_key}", fontsize=12)
    ax.legend(fontsize=10)
    plt.tight_layout()
    out = os.path.join(outdir, f"bdt_score_{signal_key}.png")
    plt.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  Saved: {out}")


def plot_roc(sig_scores, bkg_scores, signal_key, outdir):
    if not _HAS_MPL:
        return

    # Build ROC: scan threshold over the combined score range
    all_scores = np.concatenate([sig_scores, bkg_scores])
    thresholds  = np.linspace(all_scores.min(), all_scores.max(), 500)

    tpr, fpr = [], []
    for t in thresholds:
        tp = np.sum(sig_scores >= t)
        fp = np.sum(bkg_scores >= t)
        tpr.append(tp / max(len(sig_scores), 1))
        fpr.append(fp / max(len(bkg_scores), 1))

    tpr = np.array(tpr)
    fpr = np.array(fpr)
    bkg_rej = 1.0 - fpr

    # AUC via trapezoidal rule in (tpr, bkg_rej) space
    order  = np.argsort(tpr)
    auc    = np.trapz(bkg_rej[order], tpr[order])

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(tpr, bkg_rej, color="crimson", lw=2,
            label=f"ALP {signal_key}  (AUC = {auc:.4f})")
    ax.plot([0, 1], [1, 0], "k--", lw=1, label="Random")
    ax.set_xlabel("Signal efficiency", fontsize=12)
    ax.set_ylabel("Background rejection (1 − FPR)", fontsize=12)
    ax.set_title("ROC curve", fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = os.path.join(outdir, f"roc_{signal_key}.png")
    plt.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  Saved: {out}")
    print(f"  ROC AUC ({signal_key}): {auc:.4f}")


# ── Plot 5: Variable importance from TMVA XML ─────────────────────────────────

def _parse_variable_importance(xml_path):
    """
    Parse variable importance from a TMVA BDT XML weight file.

    Returns list of (variable_name, importance) sorted descending.
    """
    import xml.etree.ElementTree as ET
    tree = ET.parse(xml_path)
    root = tree.getroot()

    importances = {}
    # Variable ranking is in <Variables> or <VariableImportance>
    for vi in root.iter("VariableImportance"):
        name = vi.get("name", "")
        val  = float(vi.get("importance", 0))
        importances[name] = val

    if not importances:
        # Alternative: scan all Variable elements with a rank attribute
        for var_el in root.iter("Variable"):
            name = var_el.get("Expression", var_el.get("Label", ""))
            rank = var_el.get("VarImp", None)
            if rank is not None and name:
                importances[name] = float(rank)

    return sorted(importances.items(), key=lambda x: -x[1])


def plot_variable_importance(xml_path, signal_key, outdir):
    if not _HAS_MPL:
        return
    if not os.path.isfile(xml_path):
        print(f"  [SKIP] Weights XML not found: {xml_path}")
        return

    items = _parse_variable_importance(xml_path)
    if not items:
        print(f"  [WARN] No variable importance data found in {xml_path}")
        return

    names, vals = zip(*items)
    y_pos = np.arange(len(names))

    fig, ax = plt.subplots(figsize=(8, max(5, len(names) * 0.35)))
    bars = ax.barh(y_pos, vals, color="steelblue", edgecolor="white")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Variable importance", fontsize=11)
    ax.set_title(f"TMVA variable importance — {signal_key}", fontsize=12)

    # Annotate bar values
    for bar, v in zip(bars, vals):
        ax.text(bar.get_width() + max(vals) * 0.01, bar.get_y() + bar.get_height() / 2,
                f"{v:.4f}", va="center", fontsize=7)

    plt.tight_layout()
    out = os.path.join(outdir, f"variable_importance_{signal_key}.png")
    plt.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  Saved: {out}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Produce BDT diagnostic plots for the ALP analysis."
    )
    parser.add_argument(
        "--signal",
        required=True,
        metavar="KEY",
        help="Signal key (e.g. ma_1.0) — matches the ntuple and TMVA output file names.",
    )
    parser.add_argument(
        "--tmva-output",
        default=None,
        metavar="PATH",
        help="TMVA output ROOT file (default: weights/BDT_ALP_<signal>.root).",
    )
    parser.add_argument(
        "--method",
        default=None,
        metavar="NAME",
        help="TMVA method name used during training (default: BDT_ALP_<signal>).",
    )
    parser.add_argument(
        "--correlations-only",
        action="store_true",
        help="Only produce input variable and correlation plots (no TMVA file needed).",
    )
    parser.add_argument(
        "--nevts",
        type=int,
        default=-1,
        metavar="N",
        help="Max events to read from ntuples for plots (default: all).",
    )
    args = parser.parse_args()

    if not _HAS_MPL:
        sys.exit("ERROR: matplotlib is not available.  Install it (pip install matplotlib).")

    sig_key     = args.signal
    method_name = args.method or f"{BDT_METHOD_NAME}_{sig_key}"
    tmva_root   = args.tmva_output or os.path.join(WEIGHT_DIR, f"{method_name}.root")
    xml_path    = os.path.join(WEIGHT_DIR, f"{method_name}.weights.xml")

    os.makedirs(PLOT_DIR, exist_ok=True)

    print(f"\n{'='*55}")
    print(f"  BDT Diagnostics  —  signal: {sig_key}")
    print(f"{'='*55}\n")

    # Always produce input-variable and correlation plots
    print("Loading ntuples …")
    sig, bkg = _load_sig_bkg(sig_key, max_events=args.nevts)

    print("\n[1/5] Input variable distributions …")
    plot_input_variables(sig, bkg, sig_key, PLOT_DIR)

    print("\n[2/5] Correlation matrices …")
    plot_correlations(sig, bkg, sig_key, PLOT_DIR)

    if args.correlations_only:
        print("\nDone (correlations-only mode).")
        return

    # TMVA-output-dependent plots
    if not os.path.isfile(tmva_root):
        print(
            f"\n[SKIP] TMVA output not found: {tmva_root}\n"
            f"  Run bdt_train.py first, or pass --tmva-output <path>.\n"
            f"  (Input variable and correlation plots were still produced.)"
        )
        return

    print("\n[3/5] Reading TMVA test scores …")
    try:
        sig_scores, bkg_scores = _read_tmva_scores(tmva_root, method_name)
        print(f"  Signal test events   : {len(sig_scores)}")
        print(f"  Background test events: {len(bkg_scores)}")
    except Exception as e:
        print(f"  [ERROR] Could not read TMVA scores: {e}")
        return

    print("\n[4/5] BDT score distributions …")
    plot_bdt_score(sig_scores, bkg_scores, sig_key, method_name, PLOT_DIR)

    print("\n[5/5] ROC curve …")
    plot_roc(sig_scores, bkg_scores, sig_key, PLOT_DIR)

    print("\n[+]   Variable importance …")
    plot_variable_importance(xml_path, sig_key, PLOT_DIR)

    print(f"\n[bdt_diagnostics] All plots saved to: {PLOT_DIR}/")


if __name__ == "__main__":
    main()
