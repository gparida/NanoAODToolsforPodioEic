#!/usr/bin/env python3
"""
run_goodelectron_parallel.py

Run podio_goodElectron_analysis.py over multiple input ROOT files in parallel.
Output files are written to OUTPUT_DIR with the same filename as each input.

Usage
-----
    python3 run_goodelectron_parallel.py <INPUT_DIR> <OUTPUT_DIR> [options]

Arguments
---------
    INPUT_DIR    Directory containing hadded input ROOT files
    OUTPUT_DIR   Directory to write processed output files

Options
-------
    --pattern PATTERN   Glob pattern for input files (default: *.root)
    --jobs N            Number of files processed in parallel (default: 5)
    --nevts N           Max events to process per file
    --branchsel FILE    keep/drop branch selection file forwarded to analysis
    --dry-run           Print commands without executing

Example
-------
    python3 run_goodelectron_parallel.py \\
        /gpfs02/eic/gparidaeic/Input_PODIOfiles/hadded/10x110 \\
        /gpfs02/eic/gparidaeic/Input_PODIOfiles/processed/10x110
"""
import argparse
import glob
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

_here            = os.path.dirname(os.path.abspath(__file__))
_analysis_script = os.path.normpath(
    os.path.join(_here, "..", "python", "postprocessing", "examples",
                 "podio_goodElectron_analysis.py")
)


def run_one(input_file, output_file, nevts, branchsel, dry_run):
    """Run analysis on one file. Returns (basename, status, message)."""
    basename = os.path.basename(input_file)

    cmd = [sys.executable, _analysis_script,
           "--input",  input_file,
           "--output", output_file]
    if nevts is not None:
        cmd += ["--nevts", str(nevts)]
    if branchsel is not None:
        cmd += ["--branchsel", branchsel]

    print(f"[{basename}] Starting → {output_file}", flush=True)

    if dry_run:
        print("  [dry-run] " + " ".join(cmd))
        return basename, "ok", "dry-run"

    # Stream stdout+stderr line by line, prefixing each line with [basename]
    # so progress rates from all parallel jobs are visible in real time.
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    last_line = ""
    for line in proc.stdout:
        line = line.rstrip()
        if line:
            print(f"[{basename}] {line}", flush=True)
            last_line = line
    proc.wait()

    if proc.returncode == 0:
        size_mb = os.path.getsize(output_file) / 1e6 if os.path.exists(output_file) else 0.0
        return basename, "ok", f"{size_mb:.1f} MB written"
    else:
        return basename, "fail", last_line or f"exit code {proc.returncode}"


def main():
    parser = argparse.ArgumentParser(
        description="Run podio_goodElectron_analysis.py in parallel over many files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("input_dir",  help="Directory containing input ROOT files")
    parser.add_argument("output_dir", help="Directory for output ROOT files")
    parser.add_argument("--pattern",  default="*.root",
                        help="Glob pattern for input files (default: *.root)")
    parser.add_argument("--jobs", "-j", type=int, default=5, metavar="N",
                        help="Parallel jobs (default: 5)")
    parser.add_argument("--nevts",    type=int, default=None, metavar="N",
                        help="Max events per file")
    parser.add_argument("--branchsel", default=None, metavar="FILE",
                        help="keep/drop branch selection file")
    parser.add_argument("--dry-run",  action="store_true",
                        help="Print commands without executing")
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        sys.exit(f"ERROR: input directory not found: {args.input_dir}")
    if not os.path.isfile(_analysis_script):
        sys.exit(f"ERROR: analysis script not found: {_analysis_script}")

    input_files = sorted(glob.glob(os.path.join(args.input_dir, args.pattern)))
    if not input_files:
        sys.exit(f"No files matching '{args.pattern}' found in {args.input_dir}")

    if not args.dry_run:
        os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 60)
    print("  run_goodelectron_parallel.py")
    print("=" * 60)
    print(f"  Input dir    : {args.input_dir}")
    print(f"  Output dir   : {args.output_dir}")
    print(f"  Files found  : {len(input_files)}")
    print(f"  Parallel jobs: {args.jobs}")
    if args.nevts:
        print(f"  Max events   : {args.nevts}")
    if args.dry_run:
        print("  *** DRY RUN — no files will be written ***")
    print("=" * 60)
    for f in input_files:
        print(f"  {os.path.basename(f)}")
    print()

    jobs = [
        (f, os.path.join(args.output_dir, os.path.basename(f)))
        for f in input_files
    ]

    results = []

    if args.jobs == 1:
        for input_file, output_file in jobs:
            r = run_one(input_file, output_file, args.nevts, args.branchsel, args.dry_run)
            results.append(r)
            _, status, msg = r
            print(f"  {status.upper()}: {msg}\n")
    else:
        futures = {}
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for input_file, output_file in jobs:
                f = pool.submit(run_one, input_file, output_file,
                                args.nevts, args.branchsel, args.dry_run)
                futures[f] = os.path.basename(input_file)
            for f in as_completed(futures):
                r = f.result()
                results.append(r)
                _, status, msg = r
                print(f"  {status.upper()} [{futures[f]}]: {msg}")

    n_ok   = sum(1 for _, s, _ in results if s == "ok")
    n_fail = sum(1 for _, s, _ in results if s == "fail")

    print()
    print("=" * 60)
    print(f"  Done: {n_ok} succeeded  |  {n_fail} failed")
    print("=" * 60)

    if n_fail:
        print("\nFailed files:")
        for name, status, msg in results:
            if status == "fail":
                print(f"  {name}: {msg}")
        sys.exit(1)


if __name__ == "__main__":
    main()
