#!/usr/bin/env python3
"""
Standalone runner: NanoAODTools + PODIO file (EIC/ePIC format).

Sets up the Python path automatically so no CMSSW or installation is needed.
Requires: ROOT with PyROOT bindings (available in eic-shell).

Usage
-----
    python scripts/run_podio_example.py [OPTIONS]

Options
    --input  PATH   Path to PODIO ROOT file  (default: auto-detected)
    --nevts  N      Max events to process    (default: all)
    --first  N      Skip first N events      (default: 0)
    --list          List all collections in the file and exit
"""
import os
import sys
import argparse

# ---------------------------------------------------------------------------
# Bootstrap: register the NanoAODTools package without CMSSW.
#
# Strategy: inject namespace packages into sys.modules so that
#   from PhysicsTools.NanoAODTools.postprocessing.framework.xxx import ...
# resolves directly to NanoAODTools/python/postprocessing/framework/xxx.py
# This works regardless of whether a symlink / build directory exists.

_script_dir   = os.path.dirname(os.path.abspath(__file__))
_package_root = os.path.dirname(_script_dir)          # NanoAODTools/
_python_src   = os.path.join(_package_root, "python")  # contains postprocessing/

import types as _types

def _ensure_namespace(dotted_name, path):
    """Register a namespace package in sys.modules if not already present."""
    if dotted_name not in sys.modules:
        mod = _types.ModuleType(dotted_name)
        mod.__path__ = [path]
        mod.__package__ = dotted_name
        mod.__spec__ = None
        sys.modules[dotted_name] = mod
    return sys.modules[dotted_name]

# Build the namespace chain:
#   PhysicsTools               (virtual)
#   PhysicsTools.NanoAODTools  → NanoAODTools/python/
_phys_dir = os.path.join(_package_root, "build", "lib", "python", "PhysicsTools")
_ensure_namespace("PhysicsTools", _phys_dir)
_ensure_namespace("PhysicsTools.NanoAODTools", _python_src)

# Add python/ to sys.path so sub-imports like
#   from postprocessing.framework.podio_reader import ...
# also work when called from within the package.
if _python_src not in sys.path:
    sys.path.insert(0, _python_src)

# ---------------------------------------------------------------------------

def find_default_podio_file():
    """Walk upward looking for the PODIO file used in this project."""
    candidates = [
        os.path.join(_package_root, "..", "PODIO file", "recon_170.root"),
        os.path.join(_package_root, "..", "recon_170.root"),
        "recon_170.root",
    ]
    for c in candidates:
        p = os.path.abspath(c)
        if os.path.isfile(p):
            return p
    return None


def list_collections(filename):
    from PhysicsTools.NanoAODTools.postprocessing.framework.podio_reader import PODIOReader
    reader = PODIOReader(filename, entry_stop=1)
    print(f"\nCollections in: {filename}")
    print(f"  Total events : {reader._total_entries}")
    colls = sorted(reader.collection_names)
    print(f"  Collections  : {len(colls)}\n")
    for c in colls:
        members = reader.member_names(c)
        print(f"  {c}  [{len(members)} members]")
        for m in members[:6]:
            print(f"      .{m}")
        if len(members) > 6:
            print(f"      ... ({len(members)-6} more)")


def main():
    parser = argparse.ArgumentParser(
        description="Run NanoAODTools EIC example analysis on a PODIO file."
    )
    parser.add_argument(
        "--input", default=None,
        help="Path to PODIO ROOT file (default: auto-detected recon_170.root)"
    )
    parser.add_argument(
        "--nevts", type=int, default=None,
        help="Maximum number of events to process (default: all)"
    )
    parser.add_argument(
        "--first", type=int, default=0,
        help="Skip first N events (default: 0)"
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List all collections and exit"
    )
    args = parser.parse_args()

    # Resolve input file
    infile = args.input
    if infile is None:
        infile = find_default_podio_file()
    if infile is None:
        sys.exit(
            "ERROR: Could not find a PODIO input file.\n"
            "Pass --input /path/to/recon_170.root"
        )
    infile = os.path.abspath(infile)
    if not os.path.isfile(infile):
        sys.exit(f"ERROR: File not found: {infile}")

    print(f"Input file : {infile}")

    if args.list:
        list_collections(infile)
        return

    # Import framework (path set up above)
    from PhysicsTools.NanoAODTools.postprocessing.framework.podio_postprocessor import PODIOPostProcessor
    from PhysicsTools.NanoAODTools.postprocessing.examples.podio_example_analysis import EICExampleAnalysis

    p = PODIOPostProcessor(
        inputFiles   = [infile],
        modules      = [EICExampleAnalysis()],
        maxEntries   = args.nevts,
        firstEntry   = args.first,
        progressEvery= 100,
    )
    p.run()


if __name__ == "__main__":
    main()
