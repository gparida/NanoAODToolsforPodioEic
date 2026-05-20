#!/usr/bin/env python3
"""
PODIO equivalent of exampleAnalysis.py — run standalone without CMSSW.

Mirrors the structure of postprocessing/examples/exampleAnalysis.py:
  - Module subclass with beginJob / analyze / endJob
  - Collection access per event
  - Simple histogram (matplotlib or skipped if unavailable)
  - PostProcessor call at the bottom

Requires: ROOT with PyROOT bindings (available in eic-shell).
Run with: python3 python/postprocessing/examples/podio_exampleAnalysis.py --input /path/to/recon.root
"""
import os
import sys
import math

# ---------------------------------------------------------------------------
# Standalone path bootstrap — equivalent to `source standalone/env_standalone.sh`
# Sets up PhysicsTools.NanoAODTools namespace without any build step.
_here         = os.path.dirname(os.path.abspath(__file__))
_python_src   = os.path.dirname(os.path.dirname(_here))   # NanoAODTools/python/
_package_root = os.path.dirname(_python_src)               # NanoAODTools/

import types as _types
def _ns(name, path):
    if name not in sys.modules:
        m = _types.ModuleType(name)
        m.__path__ = [path]
        m.__package__ = name
        sys.modules[name] = m
_ns("PhysicsTools",              os.path.join(_package_root, "build", "lib", "python", "PhysicsTools"))
_ns("PhysicsTools.NanoAODTools", _python_src)
if _python_src not in sys.path:
    sys.path.insert(0, _python_src)
# ---------------------------------------------------------------------------

from PhysicsTools.NanoAODTools.postprocessing.framework.podio_eventloop     import PODIOModule
from PhysicsTools.NanoAODTools.postprocessing.framework.podio_datamodel     import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.podio_associations  import RecoMCAssociation, MCParticleNavigator
from PhysicsTools.NanoAODTools.postprocessing.framework.podio_postprocessor import PODIOPostProcessor


class PODIOExampleAnalysis(PODIOModule):
    def __init__(self):
        pass

    def beginJob(self):
        # ── Reco quantities ──────────────────────────────────────────────
        self.h_sumpt    = []   # scalar pT sum of charged particles per event

        # ── Rec↔MC association quantities ────────────────────────────────
        # NanoAOD equivalent: Jet_genJetIdx, Electron_genIdx
        self.h_match_pdg   = {}   # PDG of MC particles matched to reco tracks
        self.h_pT_res      = []   # (pT_rec - pT_mc) / pT_mc  for matched pairs

        # ── MC truth / parent navigation ─────────────────────────────────
        # NanoAOD equivalent: GenPart_genPartIdxMother, GenPart_statusFlags
        self.h_mc_stable_pdg = {}  # PDG of stable MC particles (genStatus==1)

    def analyze(self, event):
        charged  = Collection(event, "ReconstructedChargedParticles")
        neutrals = Collection(event, "ReconstructedParticles")
        mc       = Collection(event, "MCParticles")

        # ── 1. Reco-level quantities (same as before) ─────────────────────
        eventSumPx = 0.0
        eventSumPy = 0.0

        # select events with at least 1 reconstructed charged particle
        if len(charged) >= 1:
            for p in charged:
                eventSumPx += p["momentum.x"]
                eventSumPy += p["momentum.y"]
            for p in neutrals:
                eventSumPx += p["momentum.x"]
                eventSumPy += p["momentum.y"]
            self.h_sumpt.append(math.sqrt(eventSumPx**2 + eventSumPy**2))

        # ── 2. Rec↔MC truth matching via association collection ───────────
        #
        # NanoAOD:  genIdx = event.Electron_genIdx[i]   ← direct index, -1 if none
        #
        # PODIO:    assoc = RecoMCAssociation(event, "ReconstructedChargedParticleAssociations")
        #           mc_idx, weight = assoc.best_mc(i)   ← (None, 0) if unmatched
        #
        assoc = RecoMCAssociation(event, "ReconstructedChargedParticleAssociations")

        for i, rec in enumerate(charged):
            mc_idx, weight = assoc.best_mc(i)          # NanoAOD: genIdx = Electron_genIdx[i]
            if mc_idx is not None and weight > 0.5:   # quality cut on match weight
                truth = mc[mc_idx]
                pdg   = int(truth["PDG"])
                self.h_match_pdg[pdg] = self.h_match_pdg.get(pdg, 0) + 1

                # pT resolution: (reco - truth) / truth
                pT_rec = rec.pt()
                pT_mc  = math.sqrt(truth["momentum.x"]**2 + truth["momentum.y"]**2)
                if pT_mc > 0:
                    self.h_pT_res.append((pT_rec - pT_mc) / pT_mc)

        # ── 3. MC parent/daughter navigation ─────────────────────────────
        #
        # NanoAOD:  mother = GenPart[GenPart_genPartIdxMother[i]]  ← one index
        #
        # PODIO:    nav = MCParticleNavigator(event)
        #           parents = nav.parents(i)   ← list (usually 0 or 1 entries)
        #
        nav = MCParticleNavigator(event)

        for i in range(len(mc)):
            if nav.is_stable(i):                        # NanoAOD: isLastCopy flag
                pdg = int(mc[i]["PDG"])
                self.h_mc_stable_pdg[pdg] = self.h_mc_stable_pdg.get(pdg, 0) + 1

        return True

    def endJob(self):
        # ── Reco summary ─────────────────────────────────────────────────
        if self.h_sumpt:
            print(f"\nSelected {len(self.h_sumpt)} events with ≥1 charged particle.")
            print(f"  mean sumPt = {sum(self.h_sumpt)/len(self.h_sumpt):.3f} GeV")
            print(f"  max  sumPt = {max(self.h_sumpt):.3f} GeV")

        # ── Association summary ───────────────────────────────────────────
        total_matched = sum(self.h_match_pdg.values())
        print(f"\nRec↔MC matched tracks : {total_matched}")
        top_pdg = sorted(self.h_match_pdg.items(), key=lambda x: -x[1])[:8]
        print(f"  Top matched PDG codes: {top_pdg}")
        if self.h_pT_res:
            mean_res = sum(self.h_pT_res) / len(self.h_pT_res)
            print(f"  pT resolution mean   : {mean_res:.4f}  (n={len(self.h_pT_res)})")

        # ── MC truth summary ──────────────────────────────────────────────
        total_stable = sum(self.h_mc_stable_pdg.values())
        print(f"\nStable MC particles (genStatus==1): {total_stable}")
        top_stable = sorted(self.h_mc_stable_pdg.items(), key=lambda x: -x[1])[:8]
        print(f"  Top PDG codes: {top_stable}")

        # ── Histogram output ──────────────────────────────────────────────
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
            ax.hist(self.h_sumpt, bins=50, range=(0, 50))
            ax.set_xlabel(r"$\sum p_T$ [GeV]")
            ax.set_ylabel("Events")
            ax.set_title("PODIO Example Analysis")
            fig.savefig("podio_histOut.png", dpi=120)
            plt.close(fig)
            print("\nHistogram saved to podio_histOut.png")
        except ImportError:
            print("\n(matplotlib not available — skipping histogram plot)")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="PODIO example analysis — run standalone without CMSSW."
    )
    parser.add_argument("--input", default=None, help="Path to PODIO ROOT file")
    parser.add_argument("--nevts", type=int, default=None, help="Max events to process")
    args = parser.parse_args()

    # Locate input file
    infile = args.input
    if infile is None:
        _default = os.path.join(_package_root, "..", "PODIO file", "recon_170.root")
        infile = os.path.abspath(_default)
    if not os.path.isfile(infile):
        raise SystemExit("ERROR: File not found: {}\nPass --input /path/to/recon.root".format(infile))

    p = PODIOPostProcessor(
        inputFiles    = [infile],
        modules       = [PODIOExampleAnalysis()],
        maxEntries    = args.nevts,
        progressEvery = 100,
    )
    p.run()
