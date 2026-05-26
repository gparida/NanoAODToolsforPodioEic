"""
Data model classes for PODIO events (PyROOT backend).

PODIOEvent      — represents one event; entry point for collection access.
PODIOCollection — iterable collection of physics objects for one event.
PODIOObject     — single physics object with attribute-style member access.
"""
import math


class PODIOEvent:
    """
    Single-event accessor backed by a ROOT TTree entry.

    Access collections via:
        coll = event.get("ReconstructedParticles")
        n    = event.count("ReconstructedParticles")

    NanoAOD-style shorthand:
        n = event.nReconstructedParticles
    """

    def __init__(self, tree, collections, local_idx, leaf_cache=None):
        self._tree        = tree
        self._collections = collections   # {coll_name: [member, ...]}
        self._idx         = local_idx
        self._coll_cache  = {}
        self._leaf_cache  = leaf_cache or {}

    def get(self, collection_name):
        """Return a PODIOCollection for this event."""
        if collection_name not in self._coll_cache:
            if collection_name not in self._collections:
                raise KeyError("Collection '{}' not found in file.".format(collection_name))
            self._coll_cache[collection_name] = PODIOCollection(
                self._tree, collection_name, self._collections[collection_name],
                self._leaf_cache
            )
        return self._coll_cache[collection_name]

    def count(self, collection_name):
        return len(self.get(collection_name))

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        if name.startswith("n") and name[1:] in self._collections:
            return self.count(name[1:])
        raise AttributeError(
            "PODIOEvent has no attribute '{}'. "
            "Use event.get('CollectionName') to access collections.".format(name)
        )

    def __repr__(self):
        return "<PODIOEvent idx={}>".format(self._idx)


class PODIOCollection:
    """
    A collection of physics objects for one event.

    Supports len(), iteration, and index access:
        for particle in event.get("ReconstructedParticles"):
            print(particle.energy)
    """

    def __init__(self, tree, coll_name, members, leaf_cache=None):
        self._tree       = tree
        self._name       = coll_name
        self._members    = members     # list of member names
        self._data       = None        # {member: [float, ...]} — filled on first access
        self._leaf_cache = leaf_cache or {}

    def _load(self):
        """Read all member arrays for this collection from the currently loaded TTree entry."""
        if self._data is not None:
            return
        self._data = {}
        for member in self._members:
            leaf_name = "{}.{}".format(self._name, member)
            leaf = self._leaf_cache.get(leaf_name) or self._tree.GetLeaf(leaf_name)
            if leaf is not None:
                n = leaf.GetLen()
                self._data[member] = [leaf.GetValue(j) for j in range(n)]
            else:
                self._data[member] = []

    def _get_member_array(self, member):
        self._load()
        if member in self._data:
            return self._data[member]
        # Allow underscore→dot substitution: momentum_x → momentum.x
        dotted = member.replace("_", ".", 1)
        if dotted in self._data:
            return self._data[dotted]
        raise AttributeError(
            "Member '{}' not found in collection '{}'. Available: {}".format(
                member, self._name, self._members
            )
        )

    def __len__(self):
        self._load()
        if not self._data:
            return 0
        return len(next(iter(self._data.values())))

    def __getitem__(self, index):
        n = len(self)
        if index < 0:
            index += n
        if index < 0 or index >= n:
            raise IndexError(
                "Index {} out of range for '{}' (size={}).".format(index, self._name, n)
            )
        return PODIOObject(self, index)

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def array(self, member):
        """Return the per-event list of values for a member."""
        return self._get_member_array(member)

    def __repr__(self):
        return "<PODIOCollection '{}' len={}>".format(self._name, len(self))


def Collection(event, name):
    """
    Convenience wrapper — mirrors the NanoAOD Collection(event, "Muon") syntax.

    Usage:
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
        particle.momentum_x      # underscore → dot: momentum.x
        particle['momentum.x']   # explicit dot notation
        particle.p4()            # (px, py, pz, E) tuple
    """

    def __init__(self, collection, index):
        self._coll  = collection
        self._index = index

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        return self._coll._get_member_array(name)[self._index]

    def __getitem__(self, member):
        return self._coll._get_member_array(member)[self._index]

    def has(self, member):
        try:
            self._coll._get_member_array(member)
            return True
        except AttributeError:
            return False

    def p4(self):
        """Return (px, py, pz, energy) as a 4-tuple."""
        return (
            self["momentum.x"],
            self["momentum.y"],
            self["momentum.z"],
            self["energy"],
        )

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
        return "<PODIOObject '{}[{}]'>".format(self._coll._name, self._index)
