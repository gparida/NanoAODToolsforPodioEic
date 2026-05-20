"""
PODIO file reader for NanoAODTools standalone mode.

Reads EIC/ePIC PODIO ROOT files using uproot (no ROOT dependency).
PODIO collections are stored as jagged arrays under branches named:
    CollectionName/CollectionName.member
"""
import uproot
import awkward as ak


class PODIOReader:
    """
    Opens a PODIO ROOT file and provides event-by-event access.

    Usage:
        reader = PODIOReader("recon_170.root")
        for event in reader:
            particles = event.get("ReconstructedParticles")
            for p in particles:
                print(p.energy)
    """

    def __init__(self, filename, entry_start=0, entry_stop=None):
        self._file = uproot.open(filename)
        self._tree = self._file["events"]
        self._total_entries = self._tree.num_entries

        self._entry_start = entry_start
        self._entry_stop = min(
            entry_stop if entry_stop is not None else self._total_entries,
            self._total_entries,
        )

        # {coll_name: {member_name: branch_key}}
        # e.g. {"ReconstructedParticles": {"energy": "ReconstructedParticles/ReconstructedParticles.energy"}}
        self._collections = {}
        for key in self._tree.keys():
            if "/" in key:
                coll_name, branch = key.split("/", 1)
                prefix = coll_name + "."
                member = branch[len(prefix):] if branch.startswith(prefix) else branch
                self._collections.setdefault(coll_name, {})[member] = key

        # Lazy array cache: {branch_key: ak.Array}
        self._cache = {}

    # ------------------------------------------------------------------
    # Internal helpers

    def _load(self, branch_key):
        """Load (and cache) a full branch array for the selected entry range."""
        if branch_key not in self._cache:
            self._cache[branch_key] = self._tree[branch_key].array(
                entry_start=self._entry_start,
                entry_stop=self._entry_stop,
                library="ak",
            )
        return self._cache[branch_key]

    # ------------------------------------------------------------------
    # Public API

    @property
    def collection_names(self):
        return list(self._collections.keys())

    def has_collection(self, name):
        return name in self._collections

    def member_names(self, collection_name):
        return list(self._collections.get(collection_name, {}).keys())

    def __len__(self):
        return self._entry_stop - self._entry_start

    def __iter__(self):
        for local_idx in range(len(self)):
            yield PODIOEvent(self, local_idx)

    def get_event(self, local_idx):
        return PODIOEvent(self, local_idx)


# Imported here to avoid circular dependency; defined in podio_datamodel.py
from PhysicsTools.NanoAODTools.postprocessing.framework.podio_datamodel import PODIOEvent  # noqa: E402
