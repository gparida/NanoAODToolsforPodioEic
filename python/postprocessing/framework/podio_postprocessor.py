"""
High-level standalone PostProcessor for PODIO files.

Drop-in replacement for PostProcessor when working with EIC/ePIC PODIO
ROOT files. Requires PyROOT (ROOT with Python bindings), available in eic-shell.

Example
-------
    from PhysicsTools.NanoAODTools.postprocessing.framework.podio_postprocessor import PODIOPostProcessor
    from my_module import MyAnalysis

    p = PODIOPostProcessor(
        inputFiles=["recon.root"],
        modules=[MyAnalysis()],
        maxEntries=500,
    )
    p.run()
"""
import time

from PhysicsTools.NanoAODTools.postprocessing.framework.podio_reader import PODIOReader
from PhysicsTools.NanoAODTools.postprocessing.framework.podio_eventloop import podio_event_loop


class PODIOPostProcessor:
    """
    Runs a list of PODIOModule instances over one or more PODIO ROOT files.

    Parameters
    ----------
    inputFiles  : list of str — paths to PODIO ROOT files
    modules     : list of PODIOModule
    maxEntries  : int or None — maximum events to process per file
    firstEntry  : int — skip this many events at the start of each file
    progressEvery : int — print progress every N events (0 to disable)
    branchsel   : str or None — path to a keep/drop text file (same format as
                  NanoAODTools BranchSelection).  Forwarded to each module's
                  beginFile() so it can be passed to PODIOOutputWriter.
    """

    def __init__(
        self,
        inputFiles,
        modules=None,
        maxEntries=None,
        firstEntry=0,
        progressEvery=1000,
        branchsel=None,
    ):
        self.inputFiles    = inputFiles
        self.modules       = modules or []
        self.maxEntries    = maxEntries
        self.firstEntry    = firstEntry
        self.progressEvery = progressEvery
        self.branchsel     = branchsel   # keep/drop file path forwarded to beginFile

    def run(self):
        t0 = time.time()
        total_processed = 0
        total_accepted  = 0

        for m in self.modules:
            m.beginJob()

        for fname in self.inputFiles:
            print(f"\n[PODIOPostProcessor] Opening: {fname}")
            reader = PODIOReader(
                fname,
                entry_start=self.firstEntry,
                entry_stop=(
                    self.firstEntry + self.maxEntries
                    if self.maxEntries is not None
                    else None
                ),
            )
            print(
                f"  {len(reader)} events selected "
                f"(firstEntry={self.firstEntry}, maxEntries={self.maxEntries})"
            )
            print(f"  Collections available: {len(reader.collection_names)}")

            for m in self.modules:
                m.beginFile(fname, branchsel=self.branchsel,
                            intree=reader._tree)

            n_proc, n_acc, elapsed = podio_event_loop(
                self.modules,
                reader,
                max_events=self.maxEntries,
                progress_every=self.progressEvery,
            )

            for m in self.modules:
                m.endFile(fname)

            total_processed += n_proc
            total_accepted  += n_acc
            print(f"  File done in {elapsed:.1f}s.")

        for m in self.modules:
            m.endJob()

        total_time = time.time() - t0
        rate = total_processed / max(total_time, 1e-9)
        print(
            f"\n[PODIOPostProcessor] Total: {total_processed} events processed, "
            f"{total_accepted} accepted, "
            f"{total_time:.1f}s ({rate/1000:.2f} kHz)."
        )
