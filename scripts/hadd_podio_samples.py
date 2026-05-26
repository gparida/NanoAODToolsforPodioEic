#!/usr/bin/env python3
"""
hadd_podio_samples.py

hadd all recon_*.root files in each sample subdirectory into a single file.

Usage
-----
    python3 hadd_podio_samples.py <INPUT_DIR> <OUTPUT_DIR> [options]

Arguments
---------
    INPUT_DIR    Directory containing sample subdirectories (ma_0.1/, ma_0.2/, ...)
    OUTPUT_DIR   Directory under which hadded files are written, mirroring
                 the INPUT_DIR subdirectory structure

Options
-------
    --pattern PATTERN   Glob pattern for sample subdirs   (default: ma_*)
    --suffix  SUFFIX    Suffix for the output file name   (default: _hadded)
    --dry-run           Print commands without executing
    --jobs N            Run up to N hadd processes in parallel (default: 1)

Example
-------
    python3 hadd_podio_samples.py \\
        /gpfs02/eic/namjae/madgraph/out/20260417/aem_axem/out/madgraph5-3.7.0/10x110 \\
        /gpfs02/eic/gparidaeic/Input_PODIOfiles/Jae/madgraph/out/20260417/aem_axem/out/madgraph5-3.7.0/10x110

Output layout
-------------
    <OUTPUT_DIR>/ma_0.1/ma_0.1_hadded.root
    <OUTPUT_DIR>/ma_0.2/ma_0.2_hadded.root
    ...
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed


def find_samples(input_dir, pattern):
    """Return sorted list of (sample_name, [recon_files]) tuples."""
    samples = []
    for sample_dir in sorted(glob.glob(os.path.join(input_dir, pattern))):
        if not os.path.isdir(sample_dir):
            continue
        recon_files = sorted(glob.glob(os.path.join(sample_dir, "recon_*.root")))
        samples.append((os.path.basename(sample_dir), recon_files))
    return samples


def run_hadd(sample_name, recon_files, output_dir, suffix, dry_run):
    """
    hadd one sample.  Returns (sample_name, status, message).
    status: 'ok' | 'skip' | 'fail'
    """
    if not recon_files:
        return sample_name, "skip", "no recon_*.root files found"

    out_subdir = os.path.join(output_dir, sample_name)
    out_file   = os.path.join(out_subdir, f"{sample_name}{suffix}.root")

    print(f"[{sample_name}] {len(recon_files)} file(s) → {out_file}")

    if dry_run:
        print(f"  [dry-run] mkdir -p {out_subdir}")
        print(f"  [dry-run] hadd -f {out_file} " + " ".join(recon_files))
        return sample_name, "ok", "dry-run"

    os.makedirs(out_subdir, exist_ok=True)

    cmd = ["hadd", "-f", out_file] + recon_files
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        size_mb = os.path.getsize(out_file) / 1e6
        return sample_name, "ok", f"{size_mb:.1f} MB written"
    else:
        err = result.stderr.strip().splitlines()[-1] if result.stderr.strip() else "unknown error"
        return sample_name, "fail", err


def main():
    parser = argparse.ArgumentParser(
        description="hadd recon_*.root files per sample subdirectory.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("input_dir",  help="Base input directory containing sample folders")
    parser.add_argument("output_dir", help="Base output directory (mirrors input structure)")
    parser.add_argument("--pattern",  default="ma_*",
                        help="Glob pattern for sample subdirectories (default: ma_*)")
    parser.add_argument("--suffix",   default="_hadded",
                        help="Suffix appended to sample name for output file (default: _hadded)")
    parser.add_argument("--dry-run",  action="store_true",
                        help="Print commands without executing")
    parser.add_argument("--jobs", "-j", type=int, default=1, metavar="N",
                        help="Number of parallel hadd jobs (default: 1)")
    args = parser.parse_args()

    # Validate inputs
    if not os.path.isdir(args.input_dir):
        sys.exit(f"ERROR: input directory not found: {args.input_dir}")

    if not args.dry_run and shutil.which("hadd") is None:
        sys.exit("ERROR: hadd not found in PATH — load ROOT first (e.g. enter eic-shell).")

    if not args.dry_run:
        os.makedirs(args.output_dir, exist_ok=True)
        print(f"Output dir   : {args.output_dir} (created if not present)")

    samples = find_samples(args.input_dir, args.pattern)
    if not samples:
        sys.exit(f"No subdirectories matching '{args.pattern}' found in {args.input_dir}")

    print("=" * 60)
    print("  hadd_podio_samples.py")
    print("=" * 60)
    print(f"  Input dir    : {args.input_dir}")
    print(f"  Output dir   : {args.output_dir}")
    print(f"  Sample glob  : {args.pattern}")
    print(f"  Output suffix: {args.suffix}")
    print(f"  Parallel jobs: {args.jobs}")
    if args.dry_run:
        print("  *** DRY RUN — no files will be written ***")
    print("=" * 60)
    print()

    results = []

    if args.jobs == 1:
        for sample_name, recon_files in samples:
            r = run_hadd(sample_name, recon_files, args.output_dir,
                         args.suffix, args.dry_run)
            results.append(r)
            _, status, msg = r
            print(f"  {status.upper()}: {msg}\n")
    else:
        futures = {}
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for sample_name, recon_files in samples:
                f = pool.submit(run_hadd, sample_name, recon_files,
                                args.output_dir, args.suffix, args.dry_run)
                futures[f] = sample_name
            for f in as_completed(futures):
                r = f.result()
                results.append(r)
                _, status, msg = r
                print(f"  {status.upper()}: {msg}\n")

    # Summary
    n_ok   = sum(1 for _, s, _ in results if s == "ok")
    n_fail = sum(1 for _, s, _ in results if s == "fail")
    n_skip = sum(1 for _, s, _ in results if s == "skip")

    print("=" * 60)
    print(f"  Done: {n_ok} succeeded  |  {n_fail} failed  |  {n_skip} skipped")
    print("=" * 60)

    if n_fail:
        print("\nFailed samples:")
        for name, status, msg in results:
            if status == "fail":
                print(f"  {name}: {msg}")
        sys.exit(1)


if __name__ == "__main__":
    main()
