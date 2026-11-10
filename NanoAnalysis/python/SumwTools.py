"""Total genEventSumw of a set of NanoAOD files (only the Runs tree is read).

Used by SampleLooper (run.py) and by scripts/manageProduction.py, so that the number stored in
data/samples_DB.json is computed in one place only.
"""
import contextlib
import io
import math

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ZZAnalysis.NanoAnalysis.tools import get_genEventSumw


def totalSumw(files, verbose=True):
    """files: list of (chunk number, path). Returns (sum of genEventSumw, list of the chunk numbers)"""
    values = []
    for i, (n, path) in enumerate(files, 1):
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            raise IOError(f"Cannot open {path}")
        with (contextlib.nullcontext() if verbose else contextlib.redirect_stdout(io.StringIO())):
            values.append(get_genEventSumw(f))        # get_genEventSumw prints "gen= ... sumw= ..."
        f.Close()
        if verbose:
            print(f"  [{i}/{len(files)}] Chunk{n}: sumw = {values[-1]}")
    return math.fsum(values), [n for n, _ in files]
