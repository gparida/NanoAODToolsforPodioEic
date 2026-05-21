#!/usr/bin/env python3
"""
PODIO example: select good electrons and store with gen matching.

Selects reconstructed electrons passing |eta| < 4.6
and stores their kinematics plus matched gen electron kinematics.

Output file: input_base_name_processed.root with branches:
  - gElectron_pt, gElectron_eta, gElectron_phi, gElectron_E, gElectron_px, gElectron_py, gElectron_pz, gElectron_mass
  - ggenElectron_pt, ggenElectron_eta, ggenElectron_phi, ggenElectron_E, ggenElectron_px, ggenElectron_py, ggenElectron_pz, ggenElectron_mass
  - ggenElectron_genIdx (index into MCParticles, or -1 if unmatched)

Requires: ROOT with PyROOT bindings (available in eic-shell). No pip packages needed.
Run with: python3 podio_goodElectron_analysis.py --input /path/to/recon.root
"""
import os
import sys
import math

# ─────────────────────────────────────────────────────────────────────────────
# Path bootstrap (same as exampleAnalysis.py)
_here         = os.path.dirname(os.path.abspath(__file__))
_python_src   = os.path.dirname(os.path.dirname(_here))
_package_root = os.path.dirname(_python_src)

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
# ─────────────────────────────────────────────────────────────────────────────

from PhysicsTools.NanoAODTools.postprocessing.framework.podio_eventloop import PODIOModule
from PhysicsTools.NanoAODTools.postprocessing.framework.podio_datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.podio_associations import RecoMCAssociation
from PhysicsTools.NanoAODTools.postprocessing.framework.podio_postprocessor import PODIOPostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.podio_output import PODIOOutputWriter


class GoodElectronAnalysis(PODIOModule):
    """Select good electrons and store with gen matching."""

    def __init__(self):
        self.output = None
        self.n_good_electrons = 0
        self.n_gen_matched = 0

    def beginJob(self):
        print("\n[GoodElectronAnalysis] Starting...")

    def endJob(self):
        print(f"\n[GoodElectronAnalysis] Summary:")
        print(f"  Total good electrons selected: {self.n_good_electrons}")
        print(f"  Gen-matched electrons: {self.n_gen_matched}")

    def beginFile(self, input_file, output_file=None):
        """Initialize output writer and define collections."""
        if output_file is None:
            # Auto-generate output filename: input_base_processed.root
            base = os.path.splitext(input_file)[0]
            output_file = f"{base}_processed.root"

        self.output = PODIOOutputWriter(input_file, output_file)
        self.output_file = output_file

        # Define collections with kinematic variables
        variables = ["pt", "eta", "phi", "E", "px", "py", "pz", "mass"]
        self.output.define_collection("gElectron", variables)
        self.output.define_collection("ggenElectron", variables + ["genIdx"])

    def endFile(self, input_file, output_file=None):
        """Write and close output file."""
        if self.output:
            self.output.write()
            print(f"[GoodElectronAnalysis] Output written to {self.output_file}")

    def analyze(self, event):
        """Select good electrons by gen-matched PDG and fill output."""
        charged = Collection(event, "ReconstructedChargedParticles")
        mc      = Collection(event, "MCParticles")

        # Association must be loaded before the selection loop
        try:
            assoc = RecoMCAssociation(event, "ReconstructedChargedParticleAssociations")
        except (KeyError, ValueError):
            assoc = None

        # Select good electrons using gen PDG from the association (not reco PDG, which is unreliable).
        # Store (rec_idx, mc_idx) pairs so the fill loop reuses mc_idx without a second lookup.
        good_electrons = []   # list of (rec_idx, mc_idx)
        if assoc:
            for i, p in enumerate(charged):
                mc_idx, weight = assoc.best_mc(i)
                if mc_idx is None or weight <= 0.5:
                    continue
                gen_pdg = int(mc[mc_idx]["PDG"])
                if abs(gen_pdg) == 11:          # electron (11) or positron (-11)
                    if abs(p.eta()) < 4.6:
                        good_electrons.append((i, mc_idx))

        g_ele_data = {
            "pt": [], "eta": [], "phi": [], "E": [],
            "px": [], "py": [], "pz": [], "mass": [],
        }
        g_gen_ele_data = {
            "pt": [], "eta": [], "phi": [], "E": [],
            "px": [], "py": [], "pz": [], "mass": [], "genIdx": [],
        }

        for rec_idx, mc_idx in good_electrons:
            rec_e = charged[rec_idx]
            truth = mc[mc_idx]

            # Reco kinematics
            g_ele_data["pt"].append(rec_e.pt())
            g_ele_data["eta"].append(rec_e.eta())
            g_ele_data["phi"].append(rec_e.phi())
            g_ele_data["E"].append(rec_e["energy"])
            g_ele_data["px"].append(rec_e["momentum.x"])
            g_ele_data["py"].append(rec_e["momentum.y"])
            g_ele_data["pz"].append(rec_e["momentum.z"])
            g_ele_data["mass"].append(rec_e["mass"])

            # Gen kinematics — always available since mc_idx came from the selection above
            gen_px   = truth["momentum.x"]
            gen_py   = truth["momentum.y"]
            gen_pz   = truth["momentum.z"]
            gen_mass = truth["mass"]
            gen_p    = math.sqrt(gen_px**2 + gen_py**2 + gen_pz**2)
            gen_pt   = math.sqrt(gen_px**2 + gen_py**2)
            gen_E    = math.sqrt(gen_p**2 + gen_mass**2)
            gen_eta  = (0.5 * math.log((gen_p + gen_pz) / (gen_p - gen_pz))
                        if gen_p > abs(gen_pz) else math.copysign(float("inf"), gen_pz))

            g_gen_ele_data["pt"].append(gen_pt)
            g_gen_ele_data["eta"].append(gen_eta)
            g_gen_ele_data["phi"].append(math.atan2(gen_py, gen_px))
            g_gen_ele_data["E"].append(gen_E)
            g_gen_ele_data["px"].append(gen_px)
            g_gen_ele_data["py"].append(gen_py)
            g_gen_ele_data["pz"].append(gen_pz)
            g_gen_ele_data["mass"].append(truth["mass"])
            g_gen_ele_data["genIdx"].append(mc_idx)
            self.n_gen_matched += 1

        self.output.fill_collection("gElectron",    g_ele_data)
        self.output.fill_collection("ggenElectron", g_gen_ele_data)
        self.output.fill_event()
        self.n_good_electrons += len(good_electrons)

        return True


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Select good electrons via gen-matching and write output ROOT file."
    )
    parser.add_argument("--input",  default=None, help="Path to PODIO ROOT file")
    parser.add_argument("--nevts",  type=int, default=None, help="Max events to process")
    args = parser.parse_args()

    infile = args.input
    if infile is None:
        infile = os.path.abspath(os.path.join(_package_root, "..", "PODIO file", "recon_170.root"))
    if not os.path.isfile(infile):
        raise SystemExit("ERROR: File not found: {}\nPass --input /path/to/recon.root".format(infile))

    PODIOPostProcessor(
        inputFiles=[infile],
        modules=[GoodElectronAnalysis()],
        maxEntries=args.nevts,
        progressEvery=100,
    ).run()
