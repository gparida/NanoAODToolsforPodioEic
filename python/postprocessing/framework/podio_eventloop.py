"""
Module base class and event loop for PODIO standalone analysis.

PODIOModule  — base class compatible with the NanoAODTools Module interface,
               but without ROOT dependency (ROOT-free histogram output).
podio_event_loop — iterates over events, calls module hooks.
"""
import sys
import time


class PODIOModule:
    """
    Base class for PODIO analysis modules.

    Mirrors the NanoAODTools Module interface so that analysis code can
    be reused between NanoAOD and PODIO workflows with minimal changes.
    Histogramming uses standard Python (no ROOT required).

    Subclass and override beginJob / analyze / endJob:

        class MyAnalysis(PODIOModule):
            def beginJob(self):
                self.energies = []

            def analyze(self, event):
                particles = event.get("ReconstructedParticles")
                for p in particles:
                    self.energies.append(p.energy)
                return True   # keep event

            def endJob(self):
                print("Mean energy:", sum(self.energies)/len(self.energies))
    """

    def beginJob(self):
        """Called once before processing starts."""
        pass

    def endJob(self):
        """Called once after all events are processed."""
        pass

    def beginFile(self, filename, branchsel=None):
        """Called before each input file is processed.

        branchsel : str or PODIOBranchSelection or None
            Path to a keep/drop text file forwarded from PODIOPostProcessor.
            Pass it to PODIOOutputWriter(..., branchsel=branchsel) to filter
            which original branches appear in the output.
        """
        pass

    def endFile(self, filename):
        """Called after each input file is processed."""
        pass

    def analyze(self, event):
        """
        Process one event.
        Return True to accept / continue, False to reject / stop processing
        remaining modules for this event.
        """
        return True


def podio_event_loop(
    modules,
    reader,
    max_events=None,
    progress_every=1000,
    out=sys.stdout,
):
    """
    Run a list of PODIOModule instances over events in a PODIOReader.

    Parameters
    ----------
    modules      : list of PODIOModule
    reader       : PODIOReader
    max_events   : int or None — stop after this many events
    progress_every : int — print progress every N events (0 to disable)
    out          : file-like object for progress output

    Returns
    -------
    (n_processed, n_accepted, elapsed_seconds)
    """
    n_entries = len(reader)
    if max_events is not None:
        n_entries = min(n_entries, max_events)

    n_processed = 0
    n_accepted  = 0
    t0 = time.time()
    t_last = t0

    for i, event in enumerate(reader):
        if i >= n_entries:
            break

        accepted = True
        for m in modules:
            result = m.analyze(event)
            if not result:
                accepted = False
                break

        n_processed += 1
        if accepted:
            n_accepted += 1

        if progress_every and n_processed % progress_every == 0:
            t1 = time.time()
            rate = progress_every / max(t1 - t_last, 1e-9)
            out.write(
                f"  Processed {n_processed:8d}/{n_entries} events  "
                f"({100.*n_processed/n_entries:.1f}%)  "
                f"rate={rate/1000:.1f} kHz  "
                f"accepted={n_accepted}\n"
            )
            out.flush()
            t_last = t1

    elapsed = time.time() - t0
    avg_rate = n_processed / max(elapsed, 1e-9)
    out.write(
        f"Done. Processed {n_processed} events in {elapsed:.1f}s "
        f"({avg_rate/1000:.2f} kHz). Accepted {n_accepted}/{n_processed}.\n"
    )
    return n_processed, n_accepted, elapsed
