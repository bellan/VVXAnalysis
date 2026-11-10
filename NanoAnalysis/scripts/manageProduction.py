#!/usr/bin/env python3
"""
Tools to prepare the output of a sample production (tree production step), before the analysis.

In the production directory there are the chunk directories of ALL the produced samples

    EGamma12024Gv2_Chunk0/ZZ4lAnalysis.root, EGamma12024Gv2_Chunk1/..., ZZTo4l_Chunk0/, ...

plus the bookkeeping json (condorDB.json) with, for each sample, the number of jobs (= chunks) that
were submitted: the chunks of a sample are expected to be numbered 0 ... jobs-1.

Commands (see <command> --help):

  check     <dir>   only reports, for each sample, if all the chunks are present or which are missing
  organize  <dir>   moves the chunk directories of each COMPLETE sample in a directory <sample>/
  process   <dir>   as organize, but the samples whose chunks are smaller than --threshold are merged
                    in <sample>.root instead (with the function hadd of ZZAnalysis.AnalysisStep.hadd)
  transfer  <dir>   transfers the already prepared <sample>/ or <sample>.root

organize, process and transfer can also copy the result in a destination directory (--dest) and
update data/samples_DB.json (eos_path, optionally version; the old genEventSumw is removed because
it is a new production). For that, --year and --origin (MC/data) are needed to find the DB record.
A sample that is not in the DB yet is added, with empty process, color and style.

Examples:
  manageProduction.py check    PROD_2024/
  manageProduction.py organize PROD_2024/ -s ZZTo4l
  manageProduction.py process  PROD_2024/ --threshold 20 --dest /eos/user/b/bellan/samples/Nov2026 \\
                               --year 2024 --origin MC --version Nov2026 --dry-run

A sample with missing chunks is never touched: the missing chunks are listed and the exit status is 1.
"""

import argparse
import fcntl
import json
import os
import re
import shutil
import sys
import time
from dataclasses import dataclass

CHUNK_RE   = re.compile(r"^(?P<sample>.+)_Chunk(?P<n>\d+)$")
CHUNK_FILE = "ZZ4lAnalysis.root"


# ======================================================================
# Helpers
# ======================================================================
def fmtRanges(nums):
    """[0,1,2,5,7,8] -> '0-2, 5, 7-8'"""
    nums, out, i = sorted(nums), [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(str(nums[i]) if i == j else f"{nums[i]}-{nums[j]}")
        i = j + 1
    return ", ".join(out)


def dirSize(path):
    total = 0
    for root, _, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return total


def fmtSize(nbytes):
    return f"{nbytes / 1e9:.2f} GB"


def absPath(path):
    """Absolute path that keeps the symlinks the user is working through. os.path.abspath('.') uses
    os.getcwd(), that resolves them (e.g. /eos/user/b/bellan -> /eos/home-b/bellan), and ROOT/XRootD may
    not like the resolved form: here a relative path is built on $PWD, the path as the shell shows it."""
    if os.path.isabs(path):
        return os.path.normpath(path)
    pwd = os.environ.get("PWD")
    try:
        if pwd and os.path.isabs(pwd) and os.path.samefile(pwd, os.getcwd()):
            return os.path.normpath(os.path.join(pwd, path))
    except OSError:
        pass
    return os.path.abspath(path)


class inDir:
    """Context manager: run with prodDir as current directory.
    hadd.py (ZZAnalysis) works with paths RELATIVE to the production directory, which ROOT opens as
    plain local files. With absolute /eos/... paths ROOT uses the XRootD plugin (TNetXNGFile) and
    prints "[3001] Required argument not present" at every open: so the merge uses relative paths."""
    def __init__(self, path):
        self.path = path
    def __enter__(self):
        self.old = os.getcwd()
        os.chdir(self.path)
    def __exit__(self, *exc):
        os.chdir(self.old)


def remove(path):
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    else:
        os.remove(path)


# ======================================================================
# Bookkeeping and scan of the production directory
# ======================================================================
def findBookkeeping(explicit, prodDir):
    """condorDB.json: --bookkeeping if given, otherwise the directory where the script is launched,
    otherwise the production directory (the usual case is the first: the chunks are moved elsewhere)"""
    if explicit:
        return explicit
    candidates = [os.path.join(absPath("."), "condorDB.json"), os.path.join(prodDir, "condorDB.json")]
    for c in candidates:
        if os.path.isfile(c):
            print(f"bookkeeping: {c}")
            return c
    sys.exit("condorDB.json not found in " + " nor in ".join(os.path.dirname(c) for c in candidates)
             + " (use --bookkeeping)")


def loadBookkeeping(path):
    with open(path) as f:
        data = json.load(f)
    book = {}
    for e in data:
        if e["sample"] in book:
            print(f"WARNING: {e['sample']} appears more than once in {path}, the last entry is used")
        book[e["sample"]] = e
    return book


@dataclass
class Status:
    name: str
    expected: int      # number of chunks (= jobs in the bookkeeping)
    top: dict          # chunk number -> directory still in the production directory
    inside: dict       # chunk number -> directory already moved in <name>/
    hasFile: bool      # <name>.root exists
    hasDir: bool       # <name>/ exists

    def present(self):
        return set(self.top) | set(self.inside)

    def conflicts(self):
        return sorted(set(self.top) & set(self.inside))

    def missing(self):
        return sorted(set(range(self.expected)) - self.present())

    def unexpected(self):
        return sorted(n for n in self.present() if n >= self.expected)

    def badChunks(self):
        """chunks whose directory has no (or an empty) ZZ4lAnalysis.root"""
        bad = []
        for n, path in sorted({**self.inside, **self.top}.items()):
            f = os.path.join(path, CHUNK_FILE)
            if not os.path.isfile(f) or os.path.getsize(f) == 0:
                bad.append(n)
        return bad


def scan(prodDir, book):
    top = {}
    for e in os.scandir(prodDir):
        m = CHUNK_RE.match(e.name)
        if m and e.is_dir():
            top.setdefault(m["sample"], {})[int(m["n"])] = e.path

    statuses = {}
    for name, entry in book.items():
        inside = {}
        d = os.path.join(prodDir, name)
        if os.path.isdir(d):
            for sub in os.scandir(d):
                m = CHUNK_RE.match(sub.name)
                if m and m["sample"] == name and sub.is_dir():
                    inside[int(m["n"])] = sub.path
        statuses[name] = Status(name, entry["jobs"], top.get(name, {}), inside,
                                os.path.isfile(os.path.join(prodDir, f"{name}.root")),
                                os.path.isdir(d))
    unknown = sorted(set(top) - set(book))
    return statuses, unknown


def report(st):
    """Print the completeness of a sample. Returns True if the sample can be processed."""
    if st.hasFile and st.hasDir:
        print(f"[CONFLICT]   {st.name}: both {st.name}.root and the directory {st.name}/ exist, remove one of them")
        return False
    if st.hasFile:
        print(f"[MERGED]     {st.name}: {st.name}.root already exists")
        return True

    missing, bad, conf, unexp = st.missing(), st.badChunks(), st.conflicts(), st.unexpected()
    if not (missing or bad or conf or unexp):
        print(f"[OK]         {st.name}: all the {st.expected} chunks are present")
        return True

    if not st.present():
        print(f"[INCOMPLETE] {st.name}: NO chunk found, {st.expected} expected")
        return False
    print(f"[INCOMPLETE] {st.name}:")
    if missing:
        print(f"    {len(missing)} chunk(s) MISSING out of {st.expected}: {fmtRanges(missing)}")
    if bad:
        print(f"    chunk(s) without a valid {CHUNK_FILE}: {fmtRanges(bad)}")
    if conf:
        print(f"    chunk(s) both in the production directory and in {st.name}/: {fmtRanges(conf)}")
    if unexp:
        print(f"    chunk(s) not expected by the bookkeeping (>= {st.expected}): {fmtRanges(unexp)}")
    return False


# ======================================================================
# Actions on the chunks
# ======================================================================
def doOrganize(st, prodDir, dry):
    dest = os.path.join(prodDir, st.name)
    print(f"             {st.name}: moving {len(st.top)} chunk directories in {st.name}/")
    if not dry:
        os.makedirs(dest, exist_ok=True)
        for n in sorted(st.top):
            shutil.move(st.top[n], os.path.join(dest, os.path.basename(st.top[n])))
    return dest


def countEvents(path):
    import ROOT
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise IOError(f"cannot open {path}")
    t = f.Get("Events")
    if not t:
        raise IOError(f"no Events tree in {path}")
    n = t.GetEntries()
    f.Close()
    return n


def verifyMerged(inputs, merged):
    try:
        import ROOT
        ROOT.PyConfig.IgnoreCommandLineOptions = True
    except ImportError:
        print("             (ROOT not available: check of the number of events skipped)")
        return
    expected = sum(countEvents(p) for p in inputs)
    got = countEvents(merged)
    if expected != got:
        raise RuntimeError(f"merged file has {got} events, the chunks have {expected}")
    print(f"             check OK: {got} events in Events")


def doHadd(st, prodDir, dry, verify, clean):
    """Merge the chunks of a sample with the function `hadd` of ZZAnalysis.AnalysisStep.hadd (the one used by
    haddChunks.py: haddnano.py for NanoAOD, `hadd -ff` otherwise, `cp` for a single chunk), called with the
    chunk directories of THIS sample only, so the chunks are never moved:
      <sample>_Chunk*/ZZ4lAnalysis.root  ->  .<sample>.merge/ZZ4lAnalysis.root  ->  <sample>.root
    As in haddChunks.py, the paths are relative to the production directory. The output is written in a
    temporary directory and renamed only if everything is fine."""
    nums = sorted(st.top)
    finalName = f"{st.name}.root"
    final = os.path.join(prodDir, finalName)
    print(f"             {st.name}: merging {len(nums)} chunks -> {finalName}")
    if dry:
        return final
    try:
        from ZZAnalysis.AnalysisStep.hadd import hadd
    except ImportError:
        raise RuntimeError("cannot import ZZAnalysis.AnalysisStep.hadd (is the CMSSW environment set? cmsenv)")
    tmpDir = f".{st.name}.merge"
    chunkNames = [os.path.basename(st.top[n]) for n in nums]
    with inDir(prodDir):
        if os.path.exists(tmpDir):
            raise RuntimeError(f"{os.path.join(prodDir, tmpDir)} is in the way (leftover of a failed merge?)")
        merged = os.path.join(tmpDir, CHUNK_FILE)
        try:
            sys.stdout.flush()
            try:
                hadd(os.path.join(chunkNames[0], CHUNK_FILE), tmpDir, chunkNames)
            except SystemExit:                                # hadd() calls exit(1) if the command fails
                raise RuntimeError("hadd failed (see the output above)")
            if not os.path.isfile(merged) or os.path.getsize(merged) == 0:
                raise RuntimeError(f"hadd did not produce {merged}")
            if verify:
                verifyMerged([os.path.join(n, CHUNK_FILE) for n in chunkNames], merged)
            os.replace(merged, finalName)
        finally:
            if os.path.isdir(tmpDir):
                shutil.rmtree(tmpDir)
        if clean:
            os.makedirs("Chunks", exist_ok=True)
            for name in chunkNames:
                os.rename(name, os.path.join("Chunks", name))
            print(f"             chunks moved in {os.path.join(prodDir, 'Chunks')}/")
    return final


def prepare(st, args, prodDir):
    """organize / process: returns the path of <name>.root or <name>/, or None if not possible"""
    if not report(st):
        return None
    fileP, dirP = os.path.join(prodDir, f"{st.name}.root"), os.path.join(prodDir, st.name)
    if st.hasFile:
        return fileP
    if st.hasDir:                      # already (partially) organized: complete it, whatever the size
        if st.top:
            doOrganize(st, prodDir, args.dry_run)
        else:
            print(f"             {st.name}: already organized")
        return dirP

    size = sum(dirSize(p) for p in st.top.values())
    if args.cmd == "process" and size < args.threshold * 1e9:
        print(f"             {st.name}: {fmtSize(size)} < threshold ({args.threshold} GB) -> merge")
        return doHadd(st, prodDir, args.dry_run, not args.no_verify, args.clean)
    if args.cmd == "process":
        print(f"             {st.name}: {fmtSize(size)} >= threshold ({args.threshold} GB) -> chunks kept")
    return doOrganize(st, prodDir, args.dry_run)


def prepareExisting(st, prodDir):
    """transfer: the sample must be already prepared"""
    if not report(st):
        return None
    fileP, dirP = os.path.join(prodDir, f"{st.name}.root"), os.path.join(prodDir, st.name)
    if st.hasFile:
        return fileP
    if st.hasDir:
        if st.top:
            print(f"[ERROR]      {st.name}: some chunks are still outside {st.name}/, run organize first")
            return None
        return dirP
    print(f"[ERROR]      {st.name}: neither {st.name}.root nor {st.name}/ exist, run organize or process first")
    return None


# ======================================================================
# Transfer
# ======================================================================
def listing(path):
    if os.path.isfile(path):
        return {"": os.path.getsize(path)}
    out = {}
    for root, _, files in os.walk(path):
        for f in files:
            p = os.path.join(root, f)
            out[os.path.relpath(p, path)] = os.path.getsize(p)
    return out


class Progress:
    """Progress bar of a copy, proportional to the bytes copied. On a terminal the line is rewritten in place;
    otherwise (log file) a line is printed every 10%."""
    def __init__(self, label, total, nfiles):
        self.label, self.total, self.nfiles = label, max(total, 1), nfiles
        self.done, self.filesDone = 0, 0
        self.t0 = self.last = time.time()
        self.tty = sys.stdout.isatty()
        self.nextStep = 10

    def update(self, nbytes):
        self.done += nbytes
        self.show(False)

    def fileDone(self):
        self.filesDone += 1
        self.show(False)

    def show(self, final):
        now = time.time()
        pct = int(100 * self.done / self.total)
        if self.tty:
            if not final and now - self.last < 0.2:
                return
        elif not final and (pct < self.nextStep or pct >= 100):
            return
        self.last = now
        self.nextStep = (pct // 10 + 1) * 10
        elapsed = max(now - self.t0, 1e-6)
        rate = self.done / elapsed
        eta = (self.total - self.done) / rate if rate > 0 else 0
        files = f"  files {self.filesDone}/{self.nfiles}" if self.nfiles > 1 else ""
        bar = "#" * (pct // 4) + "-" * (25 - pct // 4)
        line = (f"             {self.label}: [{bar}] {pct:3d}%  {fmtSize(self.done)}/{fmtSize(self.total)}"
                f"{files}  {rate / 1e6:.0f} MB/s" + ("" if final else f"  ETA {fmtTime(eta)}"))
        if self.tty:
            print("\r" + line + "\033[K", end="\n" if final else "", flush=True)
        else:
            print(line, flush=True)

    def close(self):
        self.done = self.total
        self.show(True)


def fmtTime(sec):
    sec = int(sec)
    return f"{sec // 3600}h{sec % 3600 // 60:02d}m" if sec >= 3600 else f"{sec // 60}m{sec % 60:02d}s"


def copyFile(src, dst, prog, bufsize=16 * 1024 * 1024):
    with open(src, "rb") as fin, open(dst, "wb") as fout:
        while True:
            buf = fin.read(bufsize)
            if not buf:
                break
            fout.write(buf)
            prog.update(len(buf))
    shutil.copystat(src, dst)
    prog.fileDone()


def copyWithProgress(src, dst, label):
    """Copy a file or a directory tree, with a progress bar proportional to the size"""
    lst = listing(src)
    prog = Progress(label, sum(lst.values()), len(lst))
    if os.path.isfile(src):
        copyFile(src, dst, prog)
    else:
        for root, _, files in os.walk(src):
            target = os.path.join(dst, os.path.relpath(root, src))
            os.makedirs(target, exist_ok=True)
            for f in files:
                copyFile(os.path.join(root, f), os.path.join(target, f), prog)
    prog.close()


def transfer(product, destDir, overwrite, dry):
    base = os.path.basename(product.rstrip("/"))
    sample = base[:-5] if base.endswith(".root") else base
    final = os.path.join(destDir, base)
    tmp = os.path.join(destDir, f".{base}.partial")
    # both layouts count: <sample>.root and <sample>/ in the same place would be ambiguous for run.py
    existing = [p for p in (os.path.join(destDir, sample + ".root"), os.path.join(destDir, sample))
                if os.path.exists(p)]
    if existing and not overwrite:
        raise FileExistsError(f"already in the destination: {existing} (use --overwrite)")

    print(f"             {sample}: copying {product} -> {final}"
          + (f" (replacing {existing})" if existing else ""))
    if dry:
        return final
    os.makedirs(destDir, exist_ok=True)
    if os.path.exists(tmp):
        remove(tmp)
    copyWithProgress(product, tmp, sample)
    if listing(product) != listing(tmp):
        remove(tmp)
        raise IOError(f"copy of {product} not identical to the original (files or sizes differ)")
    for p in existing:
        remove(p)
    os.rename(tmp, final)
    return final


def removeSource(st, product, copied, prodDir, verify, dry):
    """--move: after a verified copy (and the update of the DB) remove the source: <sample>/ or <sample>.root
    and, for <sample>.root, also the chunks it was merged from (in the production directory or in Chunks/).
    Before removing the chunks the number of events of the copy in the destination is compared with the sum
    over the chunks (ROOT); if this check cannot be done or fails, the chunks are NOT removed."""
    chunks = {}
    if product.endswith(".root"):
        chunksDir = os.path.join(prodDir, "Chunks")
        if os.path.isdir(chunksDir):
            for e in os.scandir(chunksDir):
                m = CHUNK_RE.match(e.name)
                if m and m["sample"] == st.name and e.is_dir() and int(m["n"]) < st.expected:
                    chunks[int(m["n"])] = e.path
        chunks.update({n: p for n, p in st.top.items() if n < st.expected})
    if chunks and not dry and verify:
        missing = sorted(set(range(st.expected)) - set(chunks))
        if missing:
            raise RuntimeError(f"transferred, but the chunks are NOT removed: {len(missing)} chunk(s) not found "
                               f"({fmtRanges(missing)}), the number of events cannot be checked")
        try:
            import ROOT
        except ImportError:
            raise RuntimeError("transferred, but the chunks are NOT removed: ROOT not available for the check "
                               "of the number of events (use --no-verify to remove them anyway)")
        try:
            verifyMerged([os.path.join(chunks[n], CHUNK_FILE) for n in sorted(chunks)], copied)
        except Exception as e:
            raise RuntimeError(f"transferred, but the chunks are NOT removed: {e}")
    what = [product] + [chunks[n] for n in sorted(chunks)]
    if dry:
        print(f"             {st.name}: would remove the source {product}"
              + (f" and {len(chunks)} chunk directories" if chunks else ""))
        return
    for p in what:
        remove(p)
    chunksDir = os.path.join(prodDir, "Chunks")
    if os.path.isdir(chunksDir) and not os.listdir(chunksDir):
        os.rmdir(chunksDir)
    print(f"             {st.name}: removed the source {product}"
          + (f" and {len(chunks)} chunk directories" if chunks else ""))


# ======================================================================
# samples_DB.json
# ======================================================================
def findRecord(root, origin, year, name):
    for block in root["DB"]:
        if block["origin"] == origin and block["year"] == year:
            for s in block["samples"]:
                if s["name"] == name:
                    return s
    return None


def checkDB(dbPath, origin, year, names):
    """Checks that the DB can be read (now, not after hours of merging) and
    returns the samples that are not in it yet: they will be added."""
    with open(dbPath) as f:
        root = json.load(f)
    new = [n for n in names if findRecord(root, origin, year, n) is None]
    if new:
        print(f"Not in {dbPath} yet (origin={origin}, year={year}): they will be ADDED with empty "
              f"process, color and style: {new}")
    return new


def updateDB(dbPath, origin, year, name, destDir, version, dry):
    """Lock + re-read + atomic replace, like the framework does when it writes the DB
    (same lock file). NB: the file is rewritten with tab indentation.
    If the sample is not in the DB it is added; process, color and style are left empty."""
    with open(dbPath + ".lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        with open(dbPath) as f:
            root = json.load(f)
        s = findRecord(root, origin, year, name)
        changes = []
        if s is None:
            block = next((b for b in root["DB"] if b["origin"] == origin and b["year"] == year), None)
            if block is None:
                block = {"origin": origin, "year": year, "samples": []}
                root["DB"].append(block)
                changes.append(f"NEW block origin={origin}, year={year}")
            block["samples"].append({"process": "", "name": name, "eos_path": destDir,
                                     "version": version if version is not None else "",
                                     "color": None, "style": None})
            changes.append(f"NEW record: eos_path={destDir}, version={version if version is not None else ''!r}; "
                           "process, color and style left empty")
        else:
            changes.append(f"eos_path: {s.get('eos_path')} -> {destDir}")
            s["eos_path"] = destDir
            if version is not None:
                changes.append(f"version: {s.get('version')} -> {version}")
                s["version"] = version
            if "genEventSumw" in s:
                changes.append("genEventSumw removed (new production)")
                del s["genEventSumw"]
        if not dry:
            tmp = dbPath + ".tmp"
            with open(tmp, "w") as f:
                json.dump(root, f, indent="\t")
                f.write("\n")
            os.replace(tmp, dbPath)
    return changes


def defaultDB():
    base = os.environ.get("CMSSW_BASE")
    if base:
        p = os.path.join(base, "src/VVXAnalysis/NanoAnalysis/data/samples_DB.json")
        if os.path.isfile(p):
            return p
    return None


# ======================================================================
# Main
# ======================================================================
def parseArgs():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True, metavar="{check,organize,process,transfer}")

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("dir", help="production directory, with the chunk directories")
    common.add_argument("--bookkeeping", default=None, metavar="JSON",
                        help="bookkeeping json (default: ./condorDB.json, i.e. in the directory where the script is "
                             "launched, otherwise <dir>/condorDB.json)")
    common.add_argument("-s", "--sample", action="append", default=None, metavar="NAME",
                        help="only this sample (can be repeated). Default: all the samples in the bookkeeping")
    common.add_argument("--dry-run", action="store_true", help="show what would be done, change nothing")

    tr = argparse.ArgumentParser(add_help=False)
    tr.add_argument("--dest", default=None, help="copy the result (<sample>/ or <sample>.root) in this directory")
    tr.add_argument("--year", type=int, help="year of the samples (needed to update the DB)")
    tr.add_argument("--origin", choices=["MC", "data"], help="MC or data (needed to update the DB)")
    tr.add_argument("--samples-db", default=None, metavar="JSON",
                    help="samples_DB.json to update (default: the one of the CMSSW area)")
    tr.add_argument("--no-db-update", action="store_true", help="transfer only, do not touch samples_DB.json")
    tr.add_argument("--version", default=None, help="value for the 'version' field of the DB record")
    tr.add_argument("--move", action="store_true",
                    help="after the verified copy and the update of the DB, remove the source: <sample>/ or "
                         "<sample>.root and the chunks it was merged from (their number of events is checked first)")
    tr.add_argument("--no-verify", action="store_true",
                    help="do not compare the number of events: of a merged file with its chunks (process) and, "
                         "with --move, of the copy with the chunks that are going to be removed")
    tr.add_argument("--overwrite", action="store_true", help="replace the sample if already in the destination")

    sub.add_parser("check", parents=[common], help="report the missing chunks", description="Report the missing chunks")
    sub.add_parser("organize", parents=[common, tr], help="move the chunks of each complete sample in <sample>/",
                   description="Move the chunk directories of each complete sample in <sample>/")
    pr = sub.add_parser("process", parents=[common, tr],
                        help="merge the small samples in <sample>.root, organize the large ones",
                        description="hadd of the samples whose chunks are smaller than --threshold, "
                                    "<sample>/ with the chunks for the others")
    pr.add_argument("--threshold", type=float, default=10., metavar="GB",
                    help="total size of the chunks below which the sample is merged (default: %(default)s GB)")
    pr.add_argument("--clean", action="store_true", help="after a merge, move the chunks in <dir>/Chunks/")
    sub.add_parser("transfer", parents=[common, tr], help="transfer the already prepared samples",
                   description="Transfer the already prepared <sample>/ or <sample>.root")

    args = p.parse_args()
    if args.cmd == "transfer" and not args.dest:
        p.error("transfer needs --dest")
    if getattr(args, "dest", None) and not args.no_db_update:
        if not args.year or not args.origin:
            p.error("with --dest, --year and --origin are needed to update samples_DB.json "
                    "(or use --no-db-update)")
        args.samples_db = args.samples_db or defaultDB()
        if not args.samples_db:
            p.error("samples_DB.json not found: use --samples-db (or --no-db-update)")
    return args


def main():
    args = parseArgs()
    prodDir = absPath(args.dir)
    book = loadBookkeeping(findBookkeeping(args.bookkeeping, prodDir))

    names = args.sample or list(book)
    unknownReq = [n for n in names if n not in book]
    if unknownReq:
        sys.exit(f"Not in the bookkeeping: {unknownReq}")

    statuses, unknown = scan(prodDir, book)
    if unknown:
        print(f"WARNING: chunks of samples that are not in the bookkeeping are ignored: {unknown}")

    if args.cmd == "check":
        ok = [report(statuses[n]) for n in names]
        print(f"\n{sum(ok)}/{len(ok)} samples OK")
        return 0 if all(ok) else 1

    dest = absPath(args.dest) if args.dest else None
    updateDb = bool(dest) and not args.no_db_update
    if updateDb:                       # read the DB now, and list the new samples, before the long operations
        checkDB(args.samples_db, args.origin, args.year, names)

    results = {}
    for n in names:
        st = statuses[n]
        try:
            if args.cmd == "transfer":
                product = prepareExisting(st, prodDir)
            else:
                product = prepare(st, args, prodDir)
            if product is None:
                results[n] = "NOT DONE"
                continue
            results[n] = "ok"
            if dest:
                copied = transfer(product, dest, args.overwrite, args.dry_run)
                results[n] = "ok (would be transferred)" if args.dry_run else "ok, transferred"
                if updateDb:
                    for c in updateDB(args.samples_db, args.origin, args.year, n, dest, args.version, args.dry_run):
                        print(f"             DB: {c}")
                if args.move:
                    removeSource(st, product, copied, prodDir, not args.no_verify, args.dry_run)
        except Exception as e:
            print(f"[FAILED]     {n}: {e}")
            results[n] = "FAILED"

    print("\nSummary:")
    for n, r in results.items():
        print(f"  {n:<28} {r}")
    if args.dry_run:
        print("(dry run: nothing was changed)")
    return 0 if all(r.startswith("ok") for r in results.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
