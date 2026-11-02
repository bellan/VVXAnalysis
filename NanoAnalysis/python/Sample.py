import os
import re
from dataclasses import dataclass
from typing import Optional

# Name of the ROOT file inside each chunk directory:
#   <path>/<name>/<name>_Chunk<N>/ZZ4lAnalysis.root
CHUNK_FILE_NAME = "ZZ4lAnalysis.root"


@dataclass(eq=False)
class Task:
    """Unit of work: a whole-file sample (chunk=None) or one chunk of a sample"""
    sample: "Sample"
    chunk: Optional[int]
    inputPath: str

    def isChunk(self):
        return self.chunk is not None

    def label(self):
        s = self.sample
        return f"{s.year} {s.name}" + (f" Chunk{self.chunk}" if self.isChunk() else "")

    def logStem(self):
        """Path of the log, relative to logs/<analyzer>/<year>/ and without extension"""
        if self.isChunk():
            return f"{self.sample.name}/Chunk{self.chunk}"
        return self.sample.name

    def outputBase(self):
        """Base directory passed to the Histogrammers: <base>/<analyzer>_<region>/<name>.root"""
        if self.isChunk():
            return f"{self.sample.stagingDir()}/Chunk{self.chunk}"
        return self.sample.finalBase()


class Sample:
    def  __init__(self, dbblock, origin, year, local=True):

        # fixme check on origin
        
        self.year    = year
        self.origin  = origin
        self.name    = dbblock.name
        self.process = dbblock.process
        self.setLuminosity()
        if local:
            self._path = f"samples/{self.year}"
        else:
            self._path = f"{dbblock.eos_path}"
        # Cached total genEventSumw (only for samples split in chunks), as stored in the DB
        self.sumwInfo = getattr(dbblock, "genEventSumw", None)
        # Used only if both <name>.root and the directory <name>/ exist: 'file' or 'chunks'
        self.preferredLayout = None
        #fixme: to be expandend
        

    def isMC(self):
        return self.origin == 'MC'

    def __str__(self):
        return f"origin={self.origin}, year={self.year}, name={self.name}, process={self.process}"

    def __repr__(self):
        return self.__str__()

    
    def path(self):
        """Path of the single (hadd-ed) file"""
        return f"{self._path}/{self.name}.root"

    def chunkDir(self):
        return f"{self._path}/{self.name}"

    # ------------------------------------------------------------------
    # Input layout: <name>.root  (single file)  or  <name>/ (chunks)
    # ------------------------------------------------------------------
    def layout(self):
        """'file' or 'chunks'. Raises FileNotFoundError if neither exists, ValueError if both"""
        hasFile = os.path.isfile(self.path())
        hasDir  = os.path.isdir(self.chunkDir())
        if hasFile and hasDir:
            if self.preferredLayout is not None:
                return self.preferredLayout
            raise ValueError(f"{self}: both {self.path()} and the directory {self.chunkDir()} exist. "
                             "Choose with --no-split (use the single file) or --split (use the chunks), "
                             "or remove one of them")
        if hasFile:
            return "file"
        if hasDir:
            return "chunks"
        raise FileNotFoundError(f"neither {self.path()} nor the directory {self.chunkDir()} exist")

    def chunkFiles(self):
        """Sorted list of (chunk number, root file path)"""
        pattern = re.compile(rf"^{re.escape(self.name)}_Chunk(\d+)$")
        found = {}
        for d in os.listdir(self.chunkDir()):
            m = pattern.match(d)
            if m:
                found[int(m.group(1))] = os.path.join(self.chunkDir(), d, CHUNK_FILE_NAME)
        if not found:
            raise ValueError(f"{self}: no {self.name}_Chunk<N> directory in {self.chunkDir()}")
        missing = [p for p in found.values() if not os.path.isfile(p)]
        if missing:
            raise ValueError(f"{self}: {len(missing)} chunk(s) without {CHUNK_FILE_NAME}, "
                             f"e.g. {missing[0]}")
        return sorted(found.items())

    def chunkIds(self):
        return [n for n, _ in self.chunkFiles()]

    def tasks(self, chunks=None):
        """Tasks needed to analyze this sample. 'chunks' restricts to some chunk numbers"""
        if self.layout() == "file":
            if chunks:
                raise ValueError(f"{self}: --chunk requested but the sample is a single file")
            return [Task(self, None, self.path())]

        available = self.chunkFiles()
        if chunks:
            wanted = set(chunks)
            missing = sorted(wanted - {n for n, _ in available})
            if missing:
                raise ValueError(f"{self}: chunk(s) not found: {missing}")
            available = [(n, p) for n, p in available if n in wanted]
        return [Task(self, n, p) for n, p in available]

    # ------------------------------------------------------------------
    # Total genEventSumw of a sample split in chunks (cached in the DB)
    # ------------------------------------------------------------------
    def needsSumw(self):
        return self.isMC() and self.layout() == "chunks"

    def validSumw(self):
        """Cached total genEventSumw, or None if missing or computed with a different number of chunks"""
        info = self.sumwInfo
        if not info:
            return None
        nChunks = info.get("nChunks")
        if nChunks is None and "chunks" in info:      # old format: comma separated list of chunk numbers
            nChunks = len(info["chunks"].split(","))
        if nChunks != len(self.chunkFiles()):
            return None
        return info.get("value")

    # ------------------------------------------------------------------
    # Output locations
    # ------------------------------------------------------------------
    def finalBase(self):
        return f"results/{self.year}"

    def stagingDir(self):
        """Temporary outputs of the chunks, removed after the merge"""
        return f"results/{self.year}/_chunks/{self.name}"


    def setLuminosity(self):
        
        if   self.year == 2026:
            self.luminosity = 28.06
        elif self.year == 2025:
            self.luminosity = 110.63
        elif self.year == 2024:
            self.luminosity = 109.95
        elif self.year == 2023:
            self.luminosity = 28.28
        elif self.year == 2022:
            self.luminosity = 34.76
        elif self.year == 2018:
            self.luminosity = 59.56
        elif self.year == 2017:
            self.luminosity = 42.07
        elif self.year == 2016:
            self.luminosity = 36.31
        else:
            self.luminosity = -1.

    def explainYourself(self):
        return f"sample: {self.name}, process tag: {self.process}, year: {self.year}, corresponding luminosity: {self.luminosity}/fb"

    def print(self):
        print(self.explainYourself())
