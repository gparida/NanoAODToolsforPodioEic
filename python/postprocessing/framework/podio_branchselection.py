"""
Branch keep/drop selection for PODIO output files.

Mirrors the original NanoAODTools BranchSelection exactly — same text file
format, same keep/drop/keepmatch/dropmatch keywords — but documented for
PODIO collection naming conventions.

PODIO branch naming in ROOT
---------------------------
Each PODIO collection stored as a top-level TBranch whose name is the
collection name, e.g. ``MCParticles``.  Sub-leaves use dot notation:
``MCParticles.PDG``, ``MCParticles.momentum.x``, etc.
Internal PODIO relation branches start with ``_``:
``_MCParticles_parents``, ``_ReconstructedChargedParticleAssociations_rec``.

Keep/drop file syntax
---------------------
Lines are:   keep <pattern>  |  drop <pattern>
             keepmatch <regex>  |  dropmatch <regex>

- ``keep``/``drop`` use ROOT SetBranchStatus wildcard (``*`` matches anything).
  To select an entire PODIO collection including all its sub-leaves and
  relation branches, use the collection name with a trailing ``*``:

      keep MCParticles*
      keep _MCParticles*

- ``keepmatch``/``dropmatch`` use a Python regex matched against each
  top-level branch name.

Lines beginning with ``#`` and blank lines are ignored.

Usage
-----
    from PhysicsTools.NanoAODTools.postprocessing.framework.podio_branchselection import PODIOBranchSelection

    sel = PODIOBranchSelection("podio_keep_and_drop.txt")
    # Apply to a TTree *before* CloneTree(0) to filter the output:
    sel.selectBranches(intree)
    outtree = intree.CloneTree(0)
"""
import re
try:
    _Pattern = re._pattern_type
except AttributeError:
    _Pattern = re.Pattern


class PODIOBranchSelection:
    """
    Reads a keep/drop text file and applies SetBranchStatus to a TTree.

    Identical logic to NanoAODTools BranchSelection; see module docstring
    for the PODIO-specific file format.
    """

    def __init__(self, filename):
        comment = re.compile(r"#.*")
        ops = []
        with open(filename, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line or line[0] == "#":
                    continue
                line = re.sub(comment, "", line).strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) != 2:
                    print(
                        "podio_branchselection: ignoring malformed line: {!r}".format(line)
                    )
                    continue
                op, sel = parts
                if op == "keep":
                    ops.append((sel, 1))
                elif op == "drop":
                    ops.append((sel, 0))
                elif op == "keepmatch":
                    ops.append((re.compile(r"(?:{})$".format(sel)), 1))
                elif op == "dropmatch":
                    ops.append((re.compile(r"(?:{})$".format(sel)), 0))
                else:
                    print(
                        "podio_branchselection: unknown op {!r} in line: {!r}".format(
                            op, line
                        )
                    )
        self._ops = ops

    def selectBranches(self, tree):
        """
        Enable/disable branches on *tree* according to the keep/drop rules.

        Call this on the *input* tree before ``CloneTree(0)`` so that only
        the selected branches are cloned into the output tree.
        """
        tree.SetBranchStatus("*", 1)
        branch_names = [b.GetName() for b in tree.GetListOfBranches()]
        for pat, stat in self._ops:
            if isinstance(pat, _Pattern):
                for name in branch_names:
                    if re.match(pat, name):
                        tree.SetBranchStatus(name, stat)
            else:
                tree.SetBranchStatus(pat, stat)
