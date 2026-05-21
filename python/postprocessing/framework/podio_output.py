"""
Output writer for PODIO-processed ROOT files.

Creates a standalone ROOT output file with a fresh 'events' TTree.
Does not depend on input tree structure — only writes the branches
explicitly defined by the analysis module.

Usage
-----
    # In your PODIOModule.beginFile():
    self.output = PODIOOutputWriter(input_file, output_file)
    self.output.define_collection("gElectron", ["pt", "eta", "phi", "E", "px", "py", "pz", "mass"])
    self.output.define_collection("ggenElectron", ["pt", "eta", "phi", "E", "px", "py", "pz", "mass", "genIdx"])

    # In your PODIOModule.analyze():
    self.output.fill_collection("gElectron",    {"pt": [...], "eta": [...], ...})
    self.output.fill_collection("ggenElectron", {"pt": [...], ..., "genIdx": [...]})
    self.output.fill_event()

    # In your PODIOModule.endFile():
    self.output.write()

Branch naming convention
------------------------
    n<Name>         — int, number of entries in this collection for this event
    <Name>_<var>    — float array[n<Name>], or int array for variables ending in 'Idx'/'Id'/'PDG'
"""
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from array import array


class PODIOOutputWriter:
    """
    Writes user-defined collections to a new ROOT file.

    Parameters
    ----------
    input_file  : str — path to the input PODIO file (used only for naming; not read here)
    output_file : str — path to the output ROOT file to create
    """

    def __init__(self, input_file, output_file):
        self._output_path = output_file
        ROOT.gROOT.SetBatch(True)
        self._outfile = ROOT.TFile.Open(output_file, "RECREATE")
        if not self._outfile or self._outfile.IsZombie():
            raise IOError("Cannot create output file: {}".format(output_file))
        self._outfile.cd()
        self._outtree = ROOT.TTree("events", "Processed PODIO events")

        self._collections = {}    # {coll_name: [var, ...]}
        self._count_bufs  = {}    # {count_branch_name: array('i', [0])}
        self._data_bufs   = {}    # {branch_name: array}
        self._is_int      = {}    # {branch_name: bool}

    # ------------------------------------------------------------------

    def define_collection(self, name, variables):
        """
        Declare a new collection of objects.

        Variables whose name ends with 'Idx', 'Id', or equals 'PDG' are stored as int;
        all others are stored as float.

        Parameters
        ----------
        name      : str        e.g. "gElectron"
        variables : list[str]  e.g. ["pt", "eta", "phi", "E", "px", "py", "pz", "mass"]
        """
        self._collections[name] = list(variables)
        count_name = "n{}".format(name)

        # count branch (scalar int)
        self._count_bufs[count_name] = array('i', [0])
        self._outtree.Branch(count_name,
                             self._count_bufs[count_name],
                             "{}/I".format(count_name))

        # per-variable array branches
        for var in variables:
            bname = "{}_{}".format(name, var)
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
        """
        Fill all branches for collection 'name' with this event's data.

        Parameters
        ----------
        name : str
        data : dict[str, list]   e.g. {"pt": [1.2, 3.4], "eta": [-0.5, 1.1], ...}
        """
        if name not in self._collections:
            raise ValueError(
                "Collection '{}' not defined. Call define_collection first.".format(name)
            )
        vars_  = self._collections[name]
        n      = len(data[vars_[0]]) if vars_ and vars_[0] in data else 0
        count_name = "n{}".format(name)
        self._count_bufs[count_name][0] = n

        for var in vars_:
            bname  = "{}_{}".format(name, var)
            values = data.get(var, [])
            buf    = self._data_bufs[bname]

            # Grow the C-array buffer if this event is larger than pre-allocated
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
        """Commit the current event to the output tree."""
        self._outtree.Fill()

    def write(self):
        """Write the output TTree and close the file."""
        self._outfile.cd()
        self._outtree.Write("", ROOT.TObject.kOverwrite)
        self._outfile.Close()
        print("  Output written to {}".format(self._output_path))
