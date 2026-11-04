#!/usr/bin/env python3
"""
checkProd.py

Python port of checkProd.csh.

Scans production "Chunk" folders in the current directory, checks whether
each corresponding Condor job succeeded (root file present & non-empty,
exit status OK), and archives completed chunks into a "good" folder
(default AAAOK) or reports/archives failed ones into a "bad" folder
(default AAAFAIL).

Usage:
    checkProd.py <option> [goodDir] [badDir]

Options:
    help  = show this help and exit
    dry   = only print, do not archive chunks (default is to move to AAAOK/)
    quiet = do not list jobs that are still running
    mf    = move failed jobs to AAAFAIL/
    lf    = link failed jobs (then you can run cleanup.csh /
            resubmit_Condor.csh in AAAFAIL)  # FIXME: to be rechecked

[goodDir], [badDir] = folder names to use instead of AAAOK and AAAFAIL

If a condorDB.json file (list of {"sample", "files", "chunkSize", "jobs"})
is present in the working directory, a per-sample completion report is
printed at the end, based on chunk folders ("<sample>_Chunk<N>") still not
archived into the good folder.
"""

import argparse
import glob
import json
import os
import re
import shutil
import sys

CONDOR_DB_FILE = "condorDB.json"

DEBUG = False

# Map of numeric job exit codes to a human-readable description.
# Cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/JobExitCodes
#     https://twiki.cern.ch/twiki/bin/view/CMSPublic/StandardExitCodes
EXIT_STATUS_DESCRIPTIONS = {
    84: "(missing input file)",
    85: "(error reading file)",
    92: "(failed to open file)",
    134: "(Crashed)",
    137: "(killed - probably Condor problem)",
    152: "(Exceeded CPU time)",
    -152: "(Exceeded CPU time)",
    153: "(Condor crashed, see log file)",
    -153: "(Condor crashed, see log file)",
    154: "(You cancelled the job)",
    -154: "(You cancelled the job)",
}


def debug_print(*args):
    if DEBUG:
        print(*args)


def read_lines(path):
    """Return the lines of a file, or [] if it can't be read."""
    try:
        with open(path, "r", errors="replace") as f:
            return f.readlines()
    except OSError:
        return []


def grep_lines(pattern, lines):
    """Return the lines matching a (plain or regex) pattern."""
    regex = re.compile(pattern)
    return [line for line in lines if regex.search(line)]


def error_exit(message, code=1):
    print(message)
    sys.exit(code)


def find_proc_id(chunk, id_file):
    """Look up the Condor ProcId for a chunk in log/ProcIds, zero-padded to 3 digits."""
    lines = read_lines(id_file)
    matches = grep_lines(r"^" + re.escape(chunk) + r" ", lines)
    if not matches:
        error_exit(f"ERROR {chunk} not found in {id_file}")
    proc_id_raw = matches[0].split()[1]
    return f"{int(proc_id_raw):03d}"


def resolve_outpath(chunk):
    """
    Determine the folder actually containing the job output: either the
    chunk folder itself, or wherever it was transferred to on /eos.
    """
    outpath = chunk
    batch_script = os.path.join(chunk, "batchScript.sh")
    transfer_lines = grep_lines(r"^TRANSFER_DIR", read_lines(batch_script))
    transfer_line = transfer_lines[0].rstrip("\n") if transfer_lines else ""

    if transfer_line != "TRANSFER_DIR=":
        out_files = sorted(glob.glob(os.path.join(chunk, "log", "*.out")))
        if out_files:
            last_out = out_files[-1]
            transferring = grep_lines(
                r"^Transferring output", read_lines(last_out)
            )
            if transferring:
                outpath = re.sub(
                    r"^Transferring output to: ", "", transferring[-1].rstrip("\n")
                )
    return outpath


def check_root_file(outpath):
    """Return True if the ZZ4lAnalysis.root output is missing or empty."""
    filename = os.path.join(outpath, "ZZ4lAnalysis.root")
    if not os.path.exists(filename):
        debug_print("   Missing root file in", filename)
        return True
    if os.path.getsize(filename) == 0:
        debug_print("   Empty file:", filename)
        return True
    return False


def check_exit_status(outpath):
    """
    Read <outpath>/exitStatus.txt if present.
    Returns (exitStatus, fail) exactly mirroring the original tcsh logic:
      - file exists and non-empty -> exitStatus = its content, fail = True
      - file missing              -> exitStatus = 0, fail = True
      - file exists but empty     -> exitStatus = 0, fail unchanged (False)
    """
    status_file = os.path.join(outpath, "exitStatus.txt")
    if os.path.exists(status_file) and os.path.getsize(status_file) > 0:
        with open(status_file) as f:
            content = f.read().strip()
        try:
            exit_status = int(content)
        except ValueError:
            exit_status = content
        return exit_status, True
    elif not os.path.exists(status_file):
        return 0, True
    else:
        return 0, False


def check_condor_log_failure(cluster_id, proc_id, log_file):
    """
    Look in the Condor log for this job's entries and detect failure modes
    that aren't reflected in exitStatus.txt.
    Returns (exitStatus, fail) or (None, None) if nothing matched.
    """
    if not os.path.exists(log_file):
        return None, None

    log_lines = read_lines(log_file)
    job_id_lines = [line for line in log_lines if f"{cluster_id}.{proc_id}" in line]
    joined = "".join(job_id_lines)

    if re.search(r"Job removed.*time exceeded", joined):
        return -152, True
    if re.search(
        r"The job attribute PeriodicRemove expression.*evaluated to TRUE", joined
    ):
        return -153, True
    if re.search(r"Job was aborted", joined):
        return -154, True
    return None, None


def report_and_archive_failure(chunk, opt, exit_status, cluster_id, proc_id, log_file, bad_dir):
    description = ""

    if exit_status == 0:
        # Job status unclear from exitStatus.txt - inspect the Condor log directly.
        err_files = sorted(glob.glob(os.path.join(chunk, "log", "*.err")))
        if os.path.exists(log_file):
            log_lines = read_lines(log_file)
            job_id_lines = [
                line for line in log_lines if f"{cluster_id}.{proc_id}" in line
            ]
            terminated = any("Job terminated" in line for line in job_id_lines)
            if terminated:
                if err_files:
                    last_err_lines = read_lines(err_files[-1])
                    eos_problem = any(
                        "mkdir: cannot create directory '/eos" in line
                        for line in last_err_lines
                    )
                else:
                    eos_problem = False
                if eos_problem:
                    print(chunk, ": eos transfer problem, see", err_files[-1])
                else:
                    print(chunk, ": terminated, unknown failure")
            elif opt != "quiet":
                print(chunk, ": still running (or unknown failure)")
        elif opt != "quiet":  # this should be no longer possible
            print(chunk, ": still pending (or unknown failure)")
    else:
        description = EXIT_STATUS_DESCRIPTIONS.get(exit_status, "")
        print(chunk, ": failed, exit status = ", exit_status, description)

    if opt == "mf" and exit_status != 0:
        os.makedirs(bad_dir, exist_ok=True)
        shutil.move(chunk, os.path.join(bad_dir, chunk))
    elif opt == "lf" and exit_status != 0:
        os.makedirs(bad_dir, exist_ok=True)
        os.symlink(os.path.join("..", chunk), os.path.join(bad_dir, chunk))
        condor_sub_link = os.path.join(bad_dir, "condor.sub")
        if os.path.islink(condor_sub_link) or os.path.exists(condor_sub_link):
            os.remove(condor_sub_link)
        os.symlink(os.path.join("..", "condor.sub"), condor_sub_link)


def load_condor_db(path=CONDOR_DB_FILE):
    """Load the sample -> jobs map from condorDB.json, if present."""
    if not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        print(f"WARNING: could not read {path}: {e}")
        return None


def report_sample_completion(condor_db):
    """
    For each sample in condorDB.json, count how many of its chunks are still
    sitting in the working directory (i.e. not yet archived into the good
    folder) and print the completion percentage relative to "jobs".
    Assumes chunk folders are named "<sample>_Chunk<N>".
    """
    if not condor_db:
        return

    remaining_chunks = glob.glob("*Chunk*")

    print("\nCompletion status per sample:")
    for entry in condor_db:
        sample = entry.get("sample")
        jobs = entry.get("jobs", 0)
        if not sample:
            continue

        pattern = re.compile(r"^" + re.escape(sample) + r"_Chunk\d+$")
        missing = [c for c in remaining_chunks if pattern.match(c)]
        n_missing = len(missing)

        if jobs > 0:
            pct_complete = (jobs - n_missing) / jobs * 100
        else:
            pct_complete = 0.0

        print(
            f"  {sample}: {n_missing} chunk(s) still missing out of {jobs} jobs "
            f"({pct_complete:.1f}% complete)"
        )


def main():
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Scan production chunks and archive completed ones.",
    )
    parser.add_argument("opt", nargs="?", default="", help="dry | quiet | mf | lf | help")
    parser.add_argument("gooddir", nargs="?", default="AAAOK")
    parser.add_argument("baddir", nargs="?", default="AAAFAIL")
    args = parser.parse_args()

    if args.opt == "help":
        print(__doc__.strip())
        sys.exit(0)

    opt = args.opt
    good_dir = args.gooddir
    bad_dir = args.baddir

    # Chunk -> Condor ProcId map file, written by the job creation script.
    id_file = os.path.join("log", "ProcIds")
    if not os.path.exists(id_file):
        error_exit(f"ERROR: {id_file} not found")

    # Single Condor log file shared by all chunks; ClusterId comes from its name.
    log_files = sorted(glob.glob(os.path.join("log", "*.log")))
    if log_files:
        log_file = log_files[0]
        cluster_id = os.path.splitext(os.path.basename(log_file))[0]
    else:
        error_exit(f"ERROR: {os.path.join('log', '*.log')} not found")

    debug_print("ClusterId:", cluster_id)

    for chunk in sorted(glob.glob("*Chunk*")):
        fail = False

        proc_id = find_proc_id(chunk, id_file)
        debug_print(f"   {chunk}: ProcId: {proc_id}")

        outpath = resolve_outpath(chunk)

        if check_root_file(outpath):
            fail = True

        exit_status, exit_fail = check_exit_status(outpath)
        if exit_fail:
            fail = True

        if exit_status == 0:
            log_status, log_fail = check_condor_log_failure(
                cluster_id, proc_id, log_file
            )
            if log_fail:
                exit_status = log_status
                fail = True

        if not fail:
            if opt != "dry":
                os.makedirs(good_dir, exist_ok=True)
                shutil.move(chunk, os.path.join(good_dir, chunk))
        else:
            report_and_archive_failure(
                chunk, opt, exit_status, cluster_id, proc_id, log_file, bad_dir
            )

    condor_db = load_condor_db()
    report_sample_completion(condor_db)


if __name__ == "__main__":
    main()
