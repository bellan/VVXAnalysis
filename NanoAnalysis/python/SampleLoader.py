import fcntl
import json
import os
import yaml
from dataclasses import dataclass
from typing import List, Optional

from VVXAnalysis.NanoAnalysis.Sample import Sample, chunksKey

@dataclass
class DBSampleBlock:
    """Raw data from json DB"""
    name: str
    process: str
    eos_path: str
    version: str
    color: int
    style: int
    # Optional, written by the framework (never by hand): total genEventSumw of a sample
    # split in chunks, {"value": float, "chunks": "1,2,3,..."}
    genEventSumw: Optional[dict] = None


class SampleNotFoundError(Exception):
    pass


class SampleLoader:
    def __init__(self, selected_samples, db_path='data/samples_DB.json'):
        self.db_path = db_path
        with open(db_path) as f:
            self._db = json.load(f)['DB']
        self.selected_samples = selected_samples
        
    def _check_exist_in_db(self, requested: List[str], match_key: str):
        """Check if the samples exist in the DB (with no filters on year/origin)."""
        found_in_db = set()
        for block in self._db:
            for s in block['samples']:
                found_in_db.add(s[match_key])
    
        missing = set(requested) - found_in_db
        if missing:
            raise SampleNotFoundError(
                f"The samples do not exist in the DB (please check '{match_key}'): {sorted(missing)}"
            )

    def load(self, local=True) -> List[Sample]:
        years   = self.selected_samples.years
        origins = self.selected_samples.origin
        mode    = self.selected_samples.mode
        
        if mode == 'name':
            requested = self.selected_samples.by_name
            match_key = 'name'
        elif mode == 'process':
            requested = self.selected_samples.by_process
            match_key = 'process'
        elif mode == 'all':
            requested = None
            match_key = None
        else:
            raise ValueError(f"selected mode is not valid: {mode}")

        # Check in DB, separate from selection
        if requested:
            self._check_exist_in_db(requested, match_key)
            
        samples = []
        for block in self._db:
            if block['year'] not in years:
                continue
            if block['origin'] not in origins:
                continue
            for s in block['samples']:
                db_block = DBSampleBlock(**s)
                if match_key is not None and getattr(db_block, match_key) not in requested:
                    continue
                samples.append(Sample(db_block, block['origin'], block['year'],local))

        return samples

    def store_sumw(self, sample, value, chunkIds):
        """Write the total genEventSumw of a chunked sample in the DB.

        Called only by the process that prepares the analysis (run.py locally, or the
        submission on lxplus), never by the jobs. The DB is re-read under an exclusive lock
        and replaced atomically, so a concurrent run.py cannot lose or corrupt the update.
        NB: the first write reformats the whole file (indent with tabs).
        """
        info = {"value": value, "chunks": chunksKey(chunkIds)}

        with open(self.db_path + ".lock", "w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)

            with open(self.db_path) as f:
                root = json.load(f)

            found = False
            for block in root['DB']:
                if block['year'] != sample.year or block['origin'] != sample.origin:
                    continue
                for s in block['samples']:
                    if s['name'] == sample.name:
                        s['genEventSumw'] = info
                        found = True
            if not found:
                raise SampleNotFoundError(f"{sample} not found in {self.db_path}")

            tmp = self.db_path + ".tmp"
            with open(tmp, "w") as f:
                json.dump(root, f, indent="\t")
                f.write("\n")
            os.replace(tmp, self.db_path)

        sample.sumwInfo = info
        self._db = root['DB']
