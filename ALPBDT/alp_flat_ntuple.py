"""
Flat ntuple producer for the ALP BDT analysis.

Reads PODIO reconstruction files via the NanoAODTools/PODIO framework and
writes a flat ROOT TTree (one row per event) containing event-level kinematic
features needed for TMVA BDT training.

Key object reconstruction
─────────────────────────
  Scattered electron : highest-pT particle in ReconstructedParticles with
                       |charge| > 0.5, pT > PT_MIN_ELECTRON, |η| < ETA_MAX_ELECTRON.

  Photon candidates  : particles in ReconstructedParticles with |charge| < 0.5,
                       E > E_MIN_PHOTON, |η| < ETA_MAX_PHOTON.  Sorted by E.

  DIS kinematics     : taken directly from InclusiveKinematicsElectron (x, Q², y).

Missing branches (e.g. no InclusiveKinematicsElectron) are filled with the
sentinel value SENTINEL = -999.0.  Events are always written so that signal
and background TTree sizes reflect the full sample statistics.

Usage (standalone)
──────────────────
    from alp_flat_ntuple import ALPFlatNtupleProducer
    producer = ALPFlatNtupleProducer("out.root", label=1)
    # then drive via PODIOPostProcessor (see produce_ntuples.py)
"""
import math
import os
import sys
from array import array

# ── Path bootstrap: make NanoAODTools importable without CMSSW ────────────────
import types as _types

def _bootstrap_nanoaod():
    _here    = os.path.dirname(os.path.abspath(__file__))
    _nano    = os.path.dirname(_here)   # ALPBDT/ is inside NanoAODToolsforPodioEic/
    _pyroot  = os.path.join(_nano, "python")
    _physdir = os.path.join(_nano, "build", "lib", "python", "PhysicsTools")

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

_bootstrap_nanoaod()

import ROOT                                                    # noqa: E402  (eic-shell)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)

from PhysicsTools.NanoAODTools.postprocessing.framework.podio_eventloop import PODIOModule  # noqa: E402

# Sentinel for missing / undefined values in the output TTree.
SENTINEL = -999.0

# ── Branch schema ─────────────────────────────────────────────────────────────
# (branch_name, type_code)  'i' = int32, 'f' = float32
_BRANCHES = [
    # Event class label (written by producer, used as TMVA signal/bkg flag)
    ("label",         "i"),   # 1 = signal (ALP), 0 = background (QED Compton)
    ("alp_mass_label","f"),   # ALP mass hypothesis (GeV); 0 for background
    # DIS kinematics — electron method
    ("x",             "f"),   # Bjorken x
    ("Q2",            "f"),   # Momentum transfer Q² (GeV²)
    ("y",             "f"),   # Inelasticity y
    # Scattered electron
    ("e_pt",          "f"),   # pT  (GeV)
    ("e_eta",         "f"),   # pseudo-rapidity
    ("e_phi",         "f"),   # azimuth
    ("e_E",           "f"),   # energy (GeV)
    ("e_pz",          "f"),   # longitudinal momentum (GeV)
    # Photon candidates (charge ≈ 0, sorted by E)
    ("n_gamma",       "i"),   # number of photon candidates in event
    ("g1_E",          "f"),   # leading photon energy
    ("g1_pt",         "f"),
    ("g1_eta",        "f"),
    ("g1_phi",        "f"),
    ("g2_E",          "f"),   # sub-leading photon energy
    ("g2_pt",         "f"),
    ("g2_eta",        "f"),
    ("g2_phi",        "f"),
    # Di-photon system (defined only if n_gamma ≥ 2)
    ("m_gg",          "f"),   # invariant mass of leading photon pair (GeV)
    ("pt_gg",         "f"),   # transverse momentum of diphoton system
    ("eta_gg",        "f"),   # rapidity of diphoton system
    ("phi_gg",        "f"),   # azimuth of diphoton system
    ("dR_gg",         "f"),   # ΔR between the two photons
    ("dphi_gg",       "f"),   # |Δφ| between the two photons
    ("deta_gg",       "f"),   # |Δη| between the two photons
    # Event shape
    ("n_charged",     "i"),   # charged particle multiplicity
    ("n_neutral",     "i"),   # neutral particle multiplicity
    ("HT",            "f"),   # scalar HT = Σ pT (all particles, GeV)
    ("HT_gamma_frac", "f"),   # fraction of HT carried by photons
    ("MET_pt",        "f"),   # missing transverse momentum (GeV)
    ("MET_phi",       "f"),   # MET azimuth
    ("sum_E",         "f"),   # total reconstructed energy (GeV)
    # Electron–photon angular correlations
    ("dR_eg1",        "f"),   # ΔR(e, leading photon)
    ("dR_eg2",        "f"),   # ΔR(e, sub-leading photon)
    ("m_egg",         "f"),   # invariant mass of (e + leading γ + sub-leading γ) system
]


class ALPFlatNtupleProducer(PODIOModule):
    """
    Produces a flat ROOT TTree from PODIO reconstruction files.

    Parameters
    ----------
    output_file : str
        Path for the output ROOT file.
    label : int
        Event class label: 1 for signal, 0 for background.
    alp_mass : float
        ALP mass hypothesis in GeV.  Use 0.0 for background samples.
    coll_recon : str
        Name of the all-particle reconstructed collection.
    coll_kin : str
        Name of the inclusive DIS kinematics collection.
    eta_max_e, pt_min_e : float
        Acceptance cuts for the scattered-electron candidate.
    e_min_gamma, eta_max_gamma : float
        Energy and acceptance cuts for photon candidates.
    """

    def __init__(
        self,
        output_file,
        label      = 1,
        alp_mass   = 0.0,
        coll_recon = "ReconstructedParticles",
        coll_kin   = "InclusiveKinematicsElectron",
        eta_max_e      = 4.0,
        pt_min_e       = 0.5,
        e_min_gamma    = 0.2,
        eta_max_gamma  = 4.0,
    ):
        self.output_file    = output_file
        self.label          = label
        self.alp_mass       = alp_mass
        self._coll_recon    = coll_recon
        self._coll_kin      = coll_kin
        self._eta_max_e     = eta_max_e
        self._pt_min_e      = pt_min_e
        self._e_min_gamma   = e_min_gamma
        self._eta_max_gamma = eta_max_gamma

        self._tfile = None
        self._tree  = None
        self._buf   = {}
        self.n_total   = 0
        self.n_written = 0

    # ── PODIOModule interface ─────────────────────────────────────────────────

    def beginJob(self):
        os.makedirs(os.path.dirname(os.path.abspath(self.output_file)), exist_ok=True)
        self._tfile = ROOT.TFile.Open(self.output_file, "RECREATE")
        if not self._tfile or self._tfile.IsZombie():
            raise IOError(f"Cannot create output file: {self.output_file}")
        self._tree = ROOT.TTree("events", "ALP BDT flat ntuple")
        self._tree.SetDirectory(self._tfile)

        for bname, tc in _BRANCHES:
            if tc == "i":
                self._buf[bname] = array("i", [0])
                self._tree.Branch(bname, self._buf[bname], f"{bname}/I")
            else:
                self._buf[bname] = array("f", [0.0])
                self._tree.Branch(bname, self._buf[bname], f"{bname}/F")

        print(f"[ALPFlatNtuple] Output : {self.output_file}")
        print(f"[ALPFlatNtuple] label={self.label}  alp_mass={self.alp_mass} GeV")

    def endJob(self):
        self._tfile.cd()
        self._tree.Write("", ROOT.TObject.kOverwrite)
        self._tfile.Close()
        print(
            f"[ALPFlatNtuple] Written {self.n_written}/{self.n_total} events "
            f"→ {self.output_file}"
        )

    def analyze(self, event):
        self.n_total += 1
        vals = self._defaults()

        self._fill_dis_kin(event, vals)
        self._fill_particles(event, vals)
        self._write_row(vals)
        return True

    # ── Private helpers ───────────────────────────────────────────────────────

    def _defaults(self):
        """Return a dict with every branch set to its default / sentinel value."""
        d = {}
        for bname, tc in _BRANCHES:
            if tc == "i":
                d[bname] = 0
            else:
                d[bname] = SENTINEL
        d["label"]         = self.label
        d["alp_mass_label"]= self.alp_mass
        d["n_gamma"]       = 0
        d["n_charged"]     = 0
        d["n_neutral"]     = 0
        return d

    def _fill_dis_kin(self, event, vals):
        try:
            kin = event.get(self._coll_kin)
        except KeyError:
            return
        if len(kin) > 0:
            k = kin[0]
            vals["x"]  = _safe(k, "x")
            vals["Q2"] = _safe(k, "Q2")
            vals["y"]  = _safe(k, "y")

    def _fill_particles(self, event, vals):
        try:
            rp = event.get(self._coll_recon)
        except KeyError:
            return

        charged, neutral = [], []

        for p in rp:
            px = p["momentum.x"]
            py = p["momentum.y"]
            pz = p["momentum.z"]
            E  = p["energy"]
            pt = math.sqrt(px * px + py * py)
            pm = math.sqrt(px * px + py * py + pz * pz)
            eta = _eta(pm, pz)
            phi = math.atan2(py, px)
            try:
                charge = p["charge"]
            except (KeyError, AttributeError):
                charge = 0.0
            try:
                mass = p["mass"]
            except (KeyError, AttributeError):
                mass = 0.0

            rec = dict(px=px, py=py, pz=pz, E=E, pt=pt, pm=pm,
                       eta=eta, phi=phi, mass=mass)

            if abs(charge) > 0.5:
                charged.append(rec)
            else:
                neutral.append(rec)

        vals["n_charged"] = len(charged)
        vals["n_neutral"] = len(neutral)

        # ── Scattered electron ────────────────────────────────────────
        electrons = sorted(
            [p for p in charged
             if abs(p["eta"]) < self._eta_max_e and p["pt"] > self._pt_min_e],
            key=lambda p: -p["pt"],
        )
        if electrons:
            el = electrons[0]
            vals["e_pt"]  = el["pt"]
            vals["e_eta"] = el["eta"]
            vals["e_phi"] = el["phi"]
            vals["e_E"]   = el["E"]
            vals["e_pz"]  = el["pz"]

        # ── Photon candidates ─────────────────────────────────────────
        gammas = sorted(
            [p for p in neutral
             if p["E"] > self._e_min_gamma and abs(p["eta"]) < self._eta_max_gamma],
            key=lambda p: -p["E"],
        )
        vals["n_gamma"] = len(gammas)

        if len(gammas) >= 1:
            g1 = gammas[0]
            vals["g1_E"]   = g1["E"]
            vals["g1_pt"]  = g1["pt"]
            vals["g1_eta"] = g1["eta"]
            vals["g1_phi"] = g1["phi"]

        if len(gammas) >= 2:
            g2 = gammas[1]
            vals["g2_E"]   = g2["E"]
            vals["g2_pt"]  = g2["pt"]
            vals["g2_eta"] = g2["eta"]
            vals["g2_phi"] = g2["phi"]

        # ── Di-photon system ──────────────────────────────────────────
        if len(gammas) >= 2:
            g1, g2  = gammas[0], gammas[1]
            gg_px   = g1["px"] + g2["px"]
            gg_py   = g1["py"] + g2["py"]
            gg_pz   = g1["pz"] + g2["pz"]
            gg_E    = g1["E"]  + g2["E"]
            m2      = gg_E * gg_E - (gg_px**2 + gg_py**2 + gg_pz**2)
            gg_pm   = math.sqrt(gg_px**2 + gg_py**2 + gg_pz**2)

            vals["m_gg"]    = math.sqrt(max(m2, 0.0))
            vals["pt_gg"]   = math.sqrt(gg_px**2 + gg_py**2)
            vals["eta_gg"]  = _eta(gg_pm, gg_pz)
            vals["phi_gg"]  = math.atan2(gg_py, gg_px)
            vals["dR_gg"]   = _delta_r(g1["eta"], g1["phi"], g2["eta"], g2["phi"])
            vals["dphi_gg"] = _delta_phi(g1["phi"], g2["phi"])
            vals["deta_gg"] = abs(g1["eta"] - g2["eta"])

        # ── Event shape ───────────────────────────────────────────────
        HT       = sum(p["pt"] for p in charged) + sum(p["pt"] for p in neutral)
        HT_gamma = sum(p["pt"] for p in gammas)
        sum_E    = sum(p["E"]  for p in charged) + sum(p["E"]  for p in neutral)
        met_px   = -(sum(p["px"] for p in charged) + sum(p["px"] for p in neutral))
        met_py   = -(sum(p["py"] for p in charged) + sum(p["py"] for p in neutral))

        vals["HT"]            = HT
        vals["HT_gamma_frac"] = HT_gamma / HT if HT > 0 else 0.0
        vals["MET_pt"]        = math.sqrt(met_px**2 + met_py**2)
        vals["MET_phi"]       = math.atan2(met_py, met_px)
        vals["sum_E"]         = sum_E

        # ── Electron–photon correlations ──────────────────────────────
        if electrons and len(gammas) >= 1:
            el, g1 = electrons[0], gammas[0]
            vals["dR_eg1"] = _delta_r(el["eta"], el["phi"], g1["eta"], g1["phi"])

            if len(gammas) >= 2:
                g2 = gammas[1]
                vals["dR_eg2"] = _delta_r(el["eta"], el["phi"], g2["eta"], g2["phi"])

                egg_px = el["px"] + g1["px"] + g2["px"]
                egg_py = el["py"] + g1["py"] + g2["py"]
                egg_pz = el["pz"] + g1["pz"] + g2["pz"]
                egg_E  = el["E"]  + g1["E"]  + g2["E"]
                m2_egg = egg_E**2 - (egg_px**2 + egg_py**2 + egg_pz**2)
                vals["m_egg"] = math.sqrt(max(m2_egg, 0.0))

    def _write_row(self, vals):
        for bname, tc in _BRANCHES:
            v = vals.get(bname, 0 if tc == "i" else SENTINEL)
            if tc == "i":
                self._buf[bname][0] = int(v)
            else:
                self._buf[bname][0] = float(v)
        self._tree.Fill()
        self.n_written += 1


# ── Free-function kinematic helpers ──────────────────────────────────────────

def _eta(pmag, pz):
    if pmag <= abs(pz):
        return math.copysign(1e6, pz)   # very forward particle → large |η|
    return 0.5 * math.log((pmag + pz) / (pmag - pz))


def _delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    while dphi >  math.pi: dphi -= 2.0 * math.pi
    while dphi < -math.pi: dphi += 2.0 * math.pi
    return abs(dphi)


def _delta_r(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = _delta_phi(phi1, phi2)
    return math.sqrt(deta * deta + dphi * dphi)


def _safe(obj, key, default=SENTINEL):
    try:
        v = obj[key]
        return v if v == v else default   # NaN guard
    except (KeyError, IndexError, TypeError):
        return default
