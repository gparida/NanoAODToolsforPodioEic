"""
Example EIC analysis using NanoAODTools + PODIO standalone mode.

Demonstrates:
  - Accessing reconstructed particles (ReconstructedParticles)
  - Accessing MC truth particles (MCParticles)
  - Accessing inclusive DIS kinematics (InclusiveKinematicsElectron)
  - Simple histogramming without ROOT (using Python dicts + matplotlib)

Run via:
    python scripts/run_podio_example.py
"""
import math

from PhysicsTools.NanoAODTools.postprocessing.framework.podio_eventloop import PODIOModule


class EICExampleAnalysis(PODIOModule):
    """
    Fills basic distributions from an EIC PODIO reconstruction file.

    Histograms (stored as plain Python lists for matplotlib):
      - Reconstructed particle momentum magnitude
      - Reconstructed particle transverse momentum (pT)
      - Number of reconstructed particles per event
      - Inclusive DIS kinematics: x, Q2, y  (electron method)
      - MC particle PDG codes
    """

    def beginJob(self):
        # Simple accumulator lists — fill with matplotlib later
        self.h_recon_p   = []   # |p| of each ReconstructedParticle
        self.h_recon_pt  = []   # pT
        self.h_recon_e   = []   # energy
        self.h_nrecon    = []   # multiplicity per event
        self.h_x         = []   # Bjorken x (electron method)
        self.h_Q2        = []   # Q2
        self.h_y         = []   # inelasticity y
        self.h_mc_pdg    = {}   # PDG → count
        self.n_events    = 0
        print("[EICExampleAnalysis] beginJob")

    def analyze(self, event):
        self.n_events += 1

        # ---- Reconstructed particles ----
        recon = event.get("ReconstructedParticles")
        self.h_nrecon.append(len(recon))
        for p in recon:
            self.h_recon_e.append(p.energy)
            self.h_recon_pt.append(p.pt())
            self.h_recon_p.append(p.p())

        # ---- DIS kinematics (electron method) ----
        kin = event.get("InclusiveKinematicsElectron")
        if len(kin) > 0:
            k = kin[0]   # one kinematic set per event
            self.h_x.append(k["x"])
            self.h_Q2.append(k["Q2"])
            self.h_y.append(k["y"])

        # ---- MC truth particles ----
        mc = event.get("MCParticles")
        for part in mc:
            pdg = int(part["PDG"])
            self.h_mc_pdg[pdg] = self.h_mc_pdg.get(pdg, 0) + 1

        return True   # keep all events

    def endJob(self):
        print(f"\n[EICExampleAnalysis] Results after {self.n_events} events:")
        print(f"  ReconstructedParticles per event : "
              f"mean={_safe_mean(self.h_nrecon):.2f}  "
              f"min={min(self.h_nrecon) if self.h_nrecon else 0}  "
              f"max={max(self.h_nrecon) if self.h_nrecon else 0}")
        print(f"  Recon |p|  : mean={_safe_mean(self.h_recon_p):.3f} GeV")
        print(f"  Recon pT   : mean={_safe_mean(self.h_recon_pt):.3f} GeV")
        print(f"  Recon E    : mean={_safe_mean(self.h_recon_e):.3f} GeV")
        if self.h_Q2:
            print(f"  Q2 (e-meth): mean={_safe_mean(self.h_Q2):.2f} GeV^2  "
                  f"range=[{min(self.h_Q2):.2f}, {max(self.h_Q2):.2f}]")
            print(f"  x  (e-meth): mean={_safe_mean(self.h_x):.4f}")
            print(f"  y  (e-meth): mean={_safe_mean(self.h_y):.4f}")

        # Top-10 MC PDG codes
        sorted_pdg = sorted(self.h_mc_pdg.items(), key=lambda x: -x[1])[:10]
        print(f"  Top MC PDG codes: {sorted_pdg}")

        # Optionally save plots with matplotlib
        self._save_plots()

    def _save_plots(self):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("  (matplotlib not available — skipping plots)")
            return

        fig, axes = plt.subplots(2, 3, figsize=(14, 8))
        fig.suptitle("EIC PODIO Example Analysis")

        _plot_hist(axes[0, 0], self.h_recon_p,  50, r"$|p|$ [GeV]",
                   "Reconstructed particles", "Counts")
        _plot_hist(axes[0, 1], self.h_recon_pt, 50, r"$p_T$ [GeV]",
                   "Reconstructed particles", "Counts")
        _plot_hist(axes[0, 2], self.h_nrecon,   30, "N particles/event",
                   "Reconstructed multiplicity", "Events")
        _plot_hist(axes[1, 0], self.h_Q2,       50, r"$Q^2$ [GeV$^2$]",
                   "Inclusive kinematics (e-method)", "Events", log_y=True)
        _plot_hist(axes[1, 1], self.h_x,        50, "Bjorken $x$",
                   "Inclusive kinematics (e-method)", "Events", log_x=True)
        _plot_hist(axes[1, 2], self.h_y,        50, "Inelasticity $y$",
                   "Inclusive kinematics (e-method)", "Events")

        plt.tight_layout()
        outname = "podio_example_output.png"
        plt.savefig(outname, dpi=120)
        print(f"  Plots saved to {outname}")
        plt.close(fig)


# ---------------------------------------------------------------------------
# Helpers

def _safe_mean(lst):
    return sum(lst) / len(lst) if lst else float("nan")


def _plot_hist(ax, data, bins, xlabel, title, ylabel, log_x=False, log_y=False):
    if not data:
        ax.set_title(title + " (no data)")
        return
    ax.hist(data, bins=bins)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
