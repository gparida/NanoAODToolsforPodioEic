"""
Shared configuration for the ALP BDT analysis framework.

Centralises file paths, physics object selection thresholds,
BDT feature list, and TMVA hyperparameters so every script
in this directory reads from one place.
"""
import os

# ── Directory layout ──────────────────────────────────────────────────────────
# ALPBDT/ lives inside NanoAODToolsforPodioEic/ (the git repo root).
_HERE     = os.path.dirname(os.path.abspath(__file__))   # .../NanoAODToolsforPodioEic/ALPBDT/
NANO_ROOT = os.path.dirname(_HERE)                        # .../NanoAODToolsforPodioEic/
_PODIO_DIR = os.path.join(
    NANO_ROOT, "..", "..", "Input_PODIOfiles", "Poolfiles_ganesh_May26_2026"
)

NTUPLE_DIR = os.path.join(_HERE, "ntuples")
WEIGHT_DIR = os.path.join(_HERE, "weights")
PLOT_DIR   = os.path.join(_HERE, "plots")

# ── Input PODIO files ─────────────────────────────────────────────────────────
# Signal: ALP → γγ at different ALP masses (GeV)
SIGNAL_FILES = {
    "ma_0.1":  os.path.join(_PODIO_DIR, "ma_0.1_hadded.root"),
    "ma_0.2":  os.path.join(_PODIO_DIR, "ma_0.2_hadded.root"),
    "ma_0.5":  os.path.join(_PODIO_DIR, "ma_0.5_hadded.root"),
    "ma_1.0":  os.path.join(_PODIO_DIR, "ma_1.0_hadded.root"),
    "ma_2.0":  os.path.join(_PODIO_DIR, "ma_2.0_hadded.root"),
    "ma_5.0":  os.path.join(_PODIO_DIR, "ma_5.0_hadded.root"),
    "ma_10.0": os.path.join(_PODIO_DIR, "ma_10.0_hadded.root"),
    "ma_20.0": os.path.join(_PODIO_DIR, "ma_20.0_hadded.root"),
}

# Background: QED Compton (e + γ → e + γ)
BKG_FILES = {
    "qed_compton": os.path.join(_PODIO_DIR, "qed_compton_hadd.root"),
}

# ── PODIO collection names ────────────────────────────────────────────────────
# Update these if the ePIC software version changes the collection naming.
COLL_RECON_ALL      = "ReconstructedParticles"          # all reco particles
COLL_RECON_CHARGED  = "ReconstructedChargedParticles"   # charged subset
COLL_KIN_ELECTRON   = "InclusiveKinematicsElectron"     # DIS kinematics (e-method)
COLL_MC             = "MCParticles"                     # MC truth

# ── Object selection thresholds ───────────────────────────────────────────────
ETA_MAX_ELECTRON = 4.0    # |η| acceptance for scattered electron
PT_MIN_ELECTRON  = 0.5    # pT threshold for scattered electron (GeV)
E_MIN_PHOTON     = 0.2    # minimum energy for photon candidates (GeV)
ETA_MAX_PHOTON   = 4.0    # |η| acceptance for photons
CHARGE_THRESHOLD = 0.5    # |charge| > this → charged particle

# ── BDT input features ────────────────────────────────────────────────────────
# These are the branches written by alp_flat_ntuple.py and read by bdt_train.py.
# Comment out any feature to exclude it from training without touching other scripts.
FEATURES = [
    # DIS kinematics (electron method)
    "x",
    "Q2",
    "y",
    # Scattered electron
    "e_pt",
    "e_eta",
    "e_phi",
    "e_E",
    "e_pz",
    # Photon multiplicity
    "n_gamma",
    # Leading photon
    "g1_E",
    "g1_pt",
    "g1_eta",
    "g1_phi",
    # Sub-leading photon
    "g2_E",
    "g2_pt",
    "g2_eta",
    "g2_phi",
    # Di-photon system
    "m_gg",
    "pt_gg",
    "eta_gg",
    "dR_gg",
    "dphi_gg",
    "deta_gg",
    # Event shape
    "n_charged",
    "n_neutral",
    "HT",
    "HT_gamma_frac",
    "MET_pt",
    "sum_E",
    # Electron–photon correlations
    "dR_eg1",
    "dR_eg2",
    "m_egg",
]

# Spectator variables: stored in the ntuple but NOT used as BDT inputs.
# Useful for cross-checks and mass-dependent studies.
SPECTATORS = ["alp_mass_label"]

# ── TMVA BDT hyperparameters ──────────────────────────────────────────────────
BDT_METHOD_NAME = "BDT_ALP"

BDT_OPTIONS = (
    "!H:!V"
    ":NTrees=500"
    ":MaxDepth=3"
    ":MinNodeSize=2.5%"
    ":BoostType=AdaBoost"
    ":AdaBoostBeta=0.5"
    ":UseBaggedBoost=True"
    ":BaggedSampleFraction=0.5"
    ":SeparationType=GiniIndex"
    ":nCuts=20"
    ":PruneMethod=NoPruning"
)

# TMVA DataLoader split options.
# nTrain_Signal=0 / nTrain_Background=0 → TMVA chooses the split automatically.
DATALOADER_OPTIONS = (
    "nTrain_Signal=0"
    ":nTrain_Background=0"
    ":SplitMode=Random"
    ":NormMode=NumEvents"
    ":!V"
)
