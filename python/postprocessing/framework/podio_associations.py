"""
Helpers for PODIO association collections and MC particle navigation.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
COMPARISON: NanoAOD GenPart vs PODIO MCParticles + Associations
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

NanoAOD truth matching
───────────────────────
  Everything is in one flat collection.  Links are direct integer indices
  baked into the particle branches:

    GenPart_pdgId[i], GenPart_pt[i], GenPart_genPartIdxMother[i]  ← index into same array
    Jet_genJetIdx[i]    ← index into GenJet collection
    Electron_genIdx[i]  ← index into GenPart collection

PODIO truth matching
─────────────────────
  Particles and links live in separate collections.

  1. MC truth particles → MCParticles collection
     Fields: PDG, mass, momentum.x/y/z, generatorStatus, ...
     No direct parent index — instead uses begin/end ranges:
       MCParticles.parents_begin[i], .parents_end[i]
       → slice [begin:end] of _MCParticles_parents.index
       → each entry is an index back into MCParticles

  2. Rec↔MC matching → ReconstructedChargedParticleAssociations
     One row per matched pair:
       .weight               : quality of match (1.0 = full truth match)
       _..._rec.index[i]     : index into ReconstructedChargedParticles
       _..._sim.index[i]     : index into MCParticles

  Visual layout (event with 2 associations):

    ReconstructedChargedParticles        MCParticles
    ───────────────────────────────      ──────────────
    [0] e-  px=0.1 py=0.0 pz=3.2   ←─── [3] PDG=11
    [1] π+  px=0.3 py=0.1 pz=1.0        [7] PDG=211 ←──
    [2] π-  px=−0.2 py=0.0 pz=0.8  ←─── [9] PDG=-211

    ReconstructedChargedParticleAssociations:
      row 0:  rec_idx=0, sim_idx=3,  weight=1.0
      row 1:  rec_idx=2, sim_idx=9,  weight=1.0
    (rec particle [1] has no association → rec-only track)

  Key difference from NanoAOD:
    NanoAOD : jet.genJetIdx  → single index, -1 if unmatched
    PODIO   : assoc.best_mc(rec_idx) → (mc_idx, weight) or (None, 0)
              Multiple matches possible (N-to-M mapping allowed)

  generatorStatus values (MCParticles, Pythia8 convention):
    1 = stable final-state particle (equivalent to NanoAOD isLastCopy)
    2 = decayed/fragmented
    3 = hard-scatter party (incoming partons, etc.)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""


class RecoMCAssociation:
    """
    Rec ↔ MC truth matching via a PODIO association collection.

    NanoAOD equivalent: Jet_genJetIdx, Electron_genIdx

    Usage
    -----
        assoc = RecoMCAssociation(event, "ReconstructedChargedParticleAssociations")
        charged = Collection(event, "ReconstructedChargedParticles")
        mc      = Collection(event, "MCParticles")

        for i, rec_particle in enumerate(charged):
            mc_idx, weight = assoc.best_mc(i)
            if mc_idx is not None:
                truth = mc[mc_idx]
                print(f"rec[{i}] matched to MCParticle PDG={int(truth['PDG'])} weight={weight:.2f}")
    """

    def __init__(self, event, assoc_name):
        def _load(coll_name, member):
            try:
                return list(event.get(coll_name).array(member))
            except (KeyError, AttributeError):
                return []

        # Relation branch names follow the pattern _<AssocName>_rec / _<AssocName>_sim
        rec_rel = "_{}_rec".format(assoc_name)
        sim_rel = "_{}_sim".format(assoc_name)

        weights  = _load(assoc_name, "weight")
        rec_idxs = _load(rec_rel,    "index")
        sim_idxs = _load(sim_rel,    "index")

        # Default weight 1.0 if the weight branch is absent
        if not weights:
            weights = [1.0] * len(rec_idxs)

        # Build O(1) lookup maps
        self._rec_to_mc = {}   # {rec_idx: [(mc_idx, weight), ...]}
        self._mc_to_rec = {}   # {mc_idx:  [(rec_idx, weight), ...]}
        for r, s, w in zip(rec_idxs, sim_idxs, weights):
            r, s = int(r), int(s)
            self._rec_to_mc.setdefault(r, []).append((s, float(w)))
            self._mc_to_rec.setdefault(s, []).append((r, float(w)))

    # ------------------------------------------------------------------

    def mc_matches(self, rec_idx):
        """
        All MC particles matched to reconstructed particle rec_idx.
        Returns list of (mc_idx, weight).  Empty list if unmatched.
        """
        return self._rec_to_mc.get(rec_idx, [])

    def rec_matches(self, mc_idx):
        """
        All reconstructed particles matched to MC particle mc_idx.
        Returns list of (rec_idx, weight).
        """
        return self._mc_to_rec.get(mc_idx, [])

    def best_mc(self, rec_idx):
        """
        Best MC match for a reconstructed particle (highest weight).
        Returns (mc_idx, weight) or (None, 0.0) if unmatched.

        NanoAOD equivalent: Jet_genJetIdx[i]  (returns -1 when unmatched)
        """
        matches = self._rec_to_mc.get(rec_idx, [])
        if not matches:
            return None, 0.0
        return max(matches, key=lambda x: x[1])

    def best_rec(self, mc_idx):
        """
        Best reconstructed match for an MC particle (highest weight).
        Returns (rec_idx, weight) or (None, 0.0) if unmatched.
        """
        matches = self._mc_to_rec.get(mc_idx, [])
        if not matches:
            return None, 0.0
        return max(matches, key=lambda x: x[1])

    def is_matched(self, rec_idx):
        """True if this reconstructed particle has at least one MC match."""
        return rec_idx in self._rec_to_mc

    def n_associations(self):
        """Total number of rec↔MC links in this event."""
        return sum(len(v) for v in self._rec_to_mc.values())


class MCParticleNavigator:
    """
    Parent/daughter navigation for the MCParticles collection.

    NanoAOD equivalent: GenPart_genPartIdxMother, iterating over GenPart

    In NanoAOD, GenPart_genPartIdxMother[i] directly gives the index of the
    mother particle (-1 if none).  In PODIO the same information is stored
    as begin/end pointers into a separate flat index list:

        MCParticles.parents_begin[i]  ─┐
        MCParticles.parents_end[i]    ─┴→ slice of _MCParticles_parents.index
                                              → each entry is an MCParticle index

    This class wraps that two-step lookup.

    Usage
    -----
        nav = MCParticleNavigator(event)
        mc  = Collection(event, "MCParticles")

        for i in range(len(mc)):
            part = mc[i]
            print(f"PDG={int(part['PDG'])} status={int(part['generatorStatus'])}")
            for parent in nav.parents(i):
                print(f"  ← parent PDG={int(parent['PDG'])}")
            for dau in nav.daughters(i):
                print(f"  → daughter PDG={int(dau['PDG'])}")
    """

    def __init__(self, event):
        def _load(coll_name, member):
            try:
                return [int(x) for x in event.get(coll_name).array(member)]
            except (KeyError, AttributeError):
                return []

        self._pb = _load("MCParticles", "parents_begin")
        self._pe = _load("MCParticles", "parents_end")
        self._db = _load("MCParticles", "daughters_begin")
        self._de = _load("MCParticles", "daughters_end")
        self._pi = _load("_MCParticles_parents",   "index")
        self._di = _load("_MCParticles_daughters", "index")
        self._mc = event.get("MCParticles")

    # ------------------------------------------------------------------

    def parents(self, mc_idx):
        """
        List of parent PODIOObjects for MCParticle at mc_idx.
        NanoAOD equivalent: GenPart[GenPart_genPartIdxMother[mc_idx]]
        Returns [] for beam particles (no parents).
        """
        start, end = self._pb[mc_idx], self._pe[mc_idx]
        return [self._mc[self._pi[j]] for j in range(start, end)]

    def daughters(self, mc_idx):
        """
        List of daughter PODIOObjects for MCParticle at mc_idx.
        Returns [] for stable final-state particles.
        """
        start, end = self._db[mc_idx], self._de[mc_idx]
        return [self._mc[self._di[j]] for j in range(start, end)]

    def parent_indices(self, mc_idx):
        """Return raw parent indices (ints) — faster than full PODIOObjects."""
        start, end = self._pb[mc_idx], self._pe[mc_idx]
        return self._pi[start:end]

    def daughter_indices(self, mc_idx):
        """Return raw daughter indices (ints)."""
        start, end = self._db[mc_idx], self._de[mc_idx]
        return self._di[start:end]

    def is_stable(self, mc_idx):
        """
        True if generatorStatus == 1 (stable, should reach detector).
        NanoAOD equivalent: GenPart_statusFlags & isLastCopy
        """
        return int(self._mc[mc_idx]["generatorStatus"]) == 1

    def is_from_hard_scatter(self, mc_idx):
        """True if generatorStatus == 3 (hard-scatter parton level)."""
        return int(self._mc[mc_idx]["generatorStatus"]) == 3

    def has_parents(self, mc_idx):
        return self._pe[mc_idx] > self._pb[mc_idx]

    def has_daughters(self, mc_idx):
        return self._de[mc_idx] > self._db[mc_idx]

    def __len__(self):
        return len(self._mc)
