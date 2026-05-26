"""
PODIO file reader for NanoAODTools standalone mode.

Reads EIC/ePIC PODIO ROOT files using PyROOT.
Requires: ROOT with PyROOT bindings (standard in eic-shell).
"""
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True   # must be set before ROOT parses argv
ROOT.gROOT.SetBatch(True)


class PODIOReader:
    """
    Opens a PODIO ROOT file and provides event-by-event access.

    Usage:
        reader = PODIOReader("recon.root")
        for event in reader:
            particles = event.get("ReconstructedParticles")
            for p in particles:
                print(p.energy)
    """

    def __init__(self, filename, entry_start=0, entry_stop=None):
        self._tfile = ROOT.TFile.Open(filename, "READ")
        if not self._tfile or self._tfile.IsZombie():
            raise IOError("Cannot open file: {}".format(filename))

        self._tree = self._tfile.Get("events")
        if not self._tree:
            raise IOError("Tree 'events' not found in {}".format(filename))

        self._total_entries = int(self._tree.GetEntries())
        self._entry_start   = int(entry_start)
        self._entry_stop    = min(
            int(entry_stop) if entry_stop is not None else self._total_entries,
            self._total_entries,
        )

        # Build collection map and leaf cache in one pass.
        # _leaf_cache avoids repeated TTree::GetLeaf searches (linear over 700+ branches)
        # by pre-mapping leaf_name → TLeaf* at startup.
        self._collections = {}
        self._leaf_cache  = {}
        for leaf in self._tree.GetListOfLeaves():
            lname = leaf.GetName()
            self._leaf_cache[lname] = leaf
            if "." in lname:
                dot    = lname.index(".")
                coll   = lname[:dot]
                member = lname[dot + 1:]
            else:
                coll, member = lname, ""
            self._collections.setdefault(coll, [])
            if member:
                self._collections[coll].append(member)

    # ------------------------------------------------------------------
    # Public API

    @property
    def collection_names(self):
        """Physics collections (excludes internal PODIO relation branches starting with _)."""
        return [k for k in self._collections if not k.startswith("_")]

    def has_collection(self, name):
        return name in self._collections

    def member_names(self, collection_name):
        return list(self._collections.get(collection_name, []))

    def __len__(self):
        return self._entry_stop - self._entry_start

    def __iter__(self):
        for local_idx in range(len(self)):
            self._tree.GetEntry(self._entry_start + local_idx)
            yield PODIOEvent(self._tree, self._collections, local_idx,
                             self._leaf_cache)

    def get_event(self, local_idx):
        self._tree.GetEntry(self._entry_start + local_idx)
        return PODIOEvent(self._tree, self._collections, local_idx,
                          self._leaf_cache)

    def close(self):
        self._tfile.Close()


# Imported here to avoid circular dependency; defined in podio_datamodel.py
from PhysicsTools.NanoAODTools.postprocessing.framework.podio_datamodel import PODIOEvent  # noqa: E402
