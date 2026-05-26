"""
Output writer for PODIO-processed ROOT files.

Clones the input 'events' tree (preserving all original branches) and adds
new user-defined collection branches on top.  Optionally filters which
original branches are cloned via a keep/drop text file (same format as the
original NanoAODTools BranchSelection).

Usage
-----
    # In your PODIOModule.beginFile():
    self.output = PODIOOutputWriter(input_file, output_file,
                                    branchsel="podio_keep_and_drop.txt")
    self.output.define_collection("gElectron", ["pt", "eta", "phi", "E", "px", "py", "pz", "mass"])

    # In your PODIOModule.analyze():
    self.output.fill_collection("gElectron", {"pt": [...], "eta": [...], ...})
    self.output.fill_event()

    # In your PODIOModule.endFile():
    self.output.write()

Branch naming convention
------------------------
    n<Name>         — int, number of entries for this event
    <Name>_<var>    — float array[n<Name>], or int array for vars ending in 'Idx'/'Id'/'PDG'

Keep/drop file format (branchsel)
----------------------------------
    drop *                          # drop everything first
    keep MCParticles*               # keep MCParticles + all sub-leaves
    keep _MCParticles*              # keep internal PODIO relation branches
    keep ReconstructedChargedParticles*
    keep _ReconstructedChargedParticleAssociations*

    Lines starting with # and blank lines are ignored.
    Patterns use ROOT SetBranchStatus wildcard (* = anything).
    See podio_keep_and_drop.txt for a working example.
"""
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from array import array

from PhysicsTools.NanoAODTools.postprocessing.framework.podio_branchselection import PODIOBranchSelection


class PODIOOutputWriter:
    """
    Writes a new ROOT file that contains selected original input branches plus
    new user-defined collection branches.

    Parameters
    ----------
    input_file  : str — path to the input PODIO ROOT file (cloned as-is)
    output_file : str — path to the output ROOT file to create
    entry_start : int — first entry being processed (must match PODIOPostProcessor.firstEntry)
    branchsel   : str or PODIOBranchSelection or None
                  Path to a keep/drop text file (or a pre-built PODIOBranchSelection
                  object).  When provided, the branch selection is applied to the
                  input tree *before* CloneTree(0) so only the selected original
                  branches appear in the output.  New user-defined collections
                  (define_collection) are always written regardless of branchsel.
                  When None (default) all original branches are preserved.
    """

    def __init__(self, input_file, output_file, entry_start=0, branchsel=None,
                 intree=None):
        ROOT.gROOT.SetBatch(True)

        if intree is not None:
            # Fast path: reuse the already-open reader tree.
            # The reader calls GetEntry(i) before analyze(); Fill() can use those
            # buffers directly — no second GetEntry needed per event.
            self._infile      = None
            self._intree      = intree
            self._owns_infile = False
        else:
            # Standalone path: open our own file handle.
            self._infile = ROOT.TFile.Open(input_file, "READ")
            if not self._infile or self._infile.IsZombie():
                raise IOError("Cannot open input file: {}".format(input_file))
            self._intree = self._infile.Get("events")
            if not self._intree:
                raise IOError("Tree 'events' not found in {}".format(input_file))
            self._owns_infile = True

        # Apply keep/drop branch selection on the input tree before CloneTree.
        # Only applied in standalone mode — applying it to the shared reader tree
        # would prevent the analysis from reading the dropped collections.
        if branchsel is not None and self._owns_infile:
            if isinstance(branchsel, str):
                branchsel = PODIOBranchSelection(branchsel)
            branchsel.selectBranches(self._intree)
            print("  Branch selection applied ({} rules).".format(len(branchsel._ops)))

        # Output file — clone the input tree structure with 0 entries.
        # CloneTree(0) shares branch memory buffers with _intree, so
        # calling _intree.GetEntry(i) automatically populates _outtree's
        # input branches ready for Fill().
        self._outfile = ROOT.TFile.Open(output_file, "RECREATE")
        if not self._outfile or self._outfile.IsZombie():
            raise IOError("Cannot create output file: {}".format(output_file))
        self._outfile.cd()
        self._outtree = self._intree.CloneTree(0)

        self._entry       = int(entry_start)  # tracks which input entry to load next
        self._collections = {}   # {coll_name: [var, ...]}
        self._count_bufs  = {}   # {count_branch: array('i',[0])}
        self._data_bufs   = {}   # {branch_name: array}
        self._is_int      = {}   # {branch_name: bool}
        self._output_path = output_file

    # ------------------------------------------------------------------

    def define_collection(self, name, variables):
        """
        Declare a new output collection.

        Variables whose name ends with 'Idx', 'Id', or equals 'PDG' are stored
        as int; all others as float.
        """
        self._collections[name] = list(variables)
        count_name = "n{}".format(name)

        self._count_bufs[count_name] = array('i', [0])
        self._outtree.Branch(count_name,
                             self._count_bufs[count_name],
                             "{}/I".format(count_name))

        for var in variables:
            bname   = "{}_{}".format(name, var)
            use_int = var.endswith("Idx") or var.endswith("Id") or var == "PDG"
            self._is_int[bname] = use_int
            if use_int:
                self._data_bufs[bname] = array('i', [0] * 64)
                self._outtree.Branch(bname, self._data_bufs[bname],
                                     "{}[{}]/I".format(bname, count_name))
            else:
                self._data_bufs[bname] = array('f', [0.0] * 64)
                self._outtree.Branch(bname, self._data_bufs[bname],
                                     "{}[{}]/F".format(bname, count_name))

        print("  Defined collection '{}': {}".format(name, ", ".join(variables)))

    def fill_collection(self, name, data):
        """Fill all branches for collection 'name' with this event's data."""
        if name not in self._collections:
            raise ValueError(
                "Collection '{}' not defined. Call define_collection first.".format(name)
            )
        vars_      = self._collections[name]
        n          = len(data[vars_[0]]) if vars_ and vars_[0] in data else 0
        count_name = "n{}".format(name)
        self._count_bufs[count_name][0] = n

        for var in vars_:
            bname  = "{}_{}".format(name, var)
            values = data.get(var, [])
            buf    = self._data_bufs[bname]

            if n > len(buf):
                new_size = n * 2
                if self._is_int[bname]:
                    self._data_bufs[bname] = array('i', [0] * new_size)
                else:
                    self._data_bufs[bname] = array('f', [0.0] * new_size)
                self._outtree.GetBranch(bname).SetAddress(self._data_bufs[bname])
                buf = self._data_bufs[bname]

            for j, v in enumerate(values):
                buf[j] = int(v) if self._is_int[bname] else float(v)

    def fill_event(self):
        """Fill the output tree for the current event.

        In standalone mode: loads original branches from our own file handle.
        In fast mode (shared reader tree): the reader already called GetEntry,
        so we just Fill() directly from the shared buffers.
        """
        if self._owns_infile:
            self._intree.GetEntry(self._entry)
        self._entry += 1
        self._outtree.Fill()

    def write(self):
        """Write the output tree and close the output file."""
        self._outfile.cd()
        self._outtree.Write("", ROOT.TObject.kOverwrite)
        self._outfile.Close()
        if self._owns_infile:
            self._infile.Close()
        print("  Output written to {}".format(self._output_path))
