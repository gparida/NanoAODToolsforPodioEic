"""
Data model classes for PODIO events.

PODIOEvent   — represents one event; entry point for collection access.
PODIOCollection — iterable collection of physics objects for one event.
PODIOObject  — single physics object with attribute-style member access.
"""
import awkward as ak
import math


class PODIOEvent:
    """
    Single-event accessor for a PODIO file.

    Access collections via:
        coll = event.get("ReconstructedParticles")
        n    = event.count("ReconstructedParticles")   # number of objects

    NanoAOD-style shorthand also works:
        n    = event.nReconstructedParticles
    """

    def __init__(self, reader, local_idx):
        self._reader = reader
        self._idx = local_idx          # index within the loaded entry range
        self._coll_cache = {}

    def get(self, collection_name):
        """Return a PODIOCollection for this event."""
        if collection_name not in self._coll_cache:
            if not self._reader.has_collection(collection_name):
                raise KeyError(f"Collection '{collection_name}' not found in file.")
            self._coll_cache[collection_name] = PODIOCollection(
                self._reader, self._idx, collection_name
            )
        return self._coll_cache[collection_name]

    def count(self, collection_name):
        return len(self.get(collection_name))

    # NanoAOD-like: event.nReconstructedParticles
    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        if name.startswith("n") and self._reader.has_collection(name[1:]):
            return self.count(name[1:])
        raise AttributeError(
            f"PODIOEvent has no attribute '{name}'. "
            f"Use event.get('CollectionName') to access collections."
        )

    def __repr__(self):
        return f"<PODIOEvent idx={self._idx}>"


class PODIOCollection:
    """
    A collection of physics objects for one event.

    Supports len(), iteration, and index access:
        for particle in event.get("ReconstructedParticles"):
            print(particle.energy)
    """

    def __init__(self, reader, local_idx, collection_name):
        self._reader = reader
        self._idx = local_idx
        self._name = collection_name
        self._members = reader._collections[collection_name]
        self._len = None

    # ------------------------------------------------------------------
    # Internal

    def _get_member_array(self, member):
        """Return the per-event array for a member (e.g., 'energy')."""
        if member not in self._members:
            # Try underscore→dot substitution: momentum_x → momentum.x
            dotted = member.replace("_", ".", 1)
            if dotted in self._members:
                member = dotted
            else:
                available = list(self._members.keys())
                raise AttributeError(
                    f"Member '{member}' not found in collection '{self._name}'. "
                    f"Available: {available}"
                )
        arr = self._reader._load(self._members[member])
        return arr[self._idx]

    # ------------------------------------------------------------------
    # Public API

    def __len__(self):
        if self._len is None:
            if not self._members:
                self._len = 0
            else:
                first_key = next(iter(self._members.values()))
                arr = self._reader._load(first_key)
                self._len = len(arr[self._idx])
        return self._len

    def __getitem__(self, index):
        n = len(self)
        if index < 0:
            index += n
        if index < 0 or index >= n:
            raise IndexError(
                f"Index {index} out of range for collection '{self._name}' "
                f"(size={n} in this event)."
            )
        return PODIOObject(self, index)

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def array(self, member):
        """Return the full per-event array for a member as a list."""
        return list(self._get_member_array(member))

    def __repr__(self):
        return f"<PODIOCollection '{self._name}' len={len(self)}>"


def Collection(event, name):
    """
    Convenience wrapper — mirrors the NanoAOD Collection(event, "Muon") syntax.

    Usage (identical pattern to exampleAnalysis.py):
        charged = Collection(event, "ReconstructedChargedParticles")
        for p in charged:
            print(p.energy)
    """
    return event.get(name)


class PODIOObject:
    """
    Single physics object inside a PODIOCollection.

    Member access:
        particle.energy          # direct member
        particle.momentum_x      # dot-separated: momentum.x
        particle['momentum.x']   # explicit dot notation
        particle.p4()            # TLorentzVector-like (x,y,z,E)
    """

    def __init__(self, collection, index):
        self._coll = collection
        self._index = index

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        arr = self._coll._get_member_array(name)
        val = arr[self._index]
        # Convert awkward scalar to Python native type
        try:
            return float(val)
        except (TypeError, ValueError):
            return val

    def __getitem__(self, member):
        arr = self._coll._get_member_array(member)
        val = arr[self._index]
        try:
            return float(val)
        except (TypeError, ValueError):
            return val

    def has(self, member):
        """Check if member exists in this collection."""
        members = self._coll._members
        return member in members or member.replace("_", ".", 1) in members

    def p4(self):
        """
        Return (px, py, pz, energy) as a simple 4-tuple.
        Works for collections that have momentum.x/y/z and energy members.
        """
        px = self["momentum.x"]
        py = self["momentum.y"]
        pz = self["momentum.z"]
        e  = self["energy"]
        return (px, py, pz, e)

    def pt(self):
        px = self["momentum.x"]
        py = self["momentum.y"]
        return math.sqrt(px * px + py * py)

    def p(self):
        px = self["momentum.x"]
        py = self["momentum.y"]
        pz = self["momentum.z"]
        return math.sqrt(px * px + py * py + pz * pz)

    def eta(self):
        pmag = self.p()
        pz   = self["momentum.z"]
        if pmag <= abs(pz):
            return math.copysign(float("inf"), pz)
        return 0.5 * math.log((pmag + pz) / (pmag - pz))

    def phi(self):
        return math.atan2(self["momentum.y"], self["momentum.x"])

    def __repr__(self):
        return (
            f"<PODIOObject '{self._coll._name}[{self._index}]'>"
        )
