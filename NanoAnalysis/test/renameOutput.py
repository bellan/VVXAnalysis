#! /usr/bin/env python3

from pathlib import Path
import re
import sys

inputDir  = Path(sys.argv[1])
outputDir = Path(sys.argv[2])

for directory in inputDir.iterdir():

    if not directory.is_dir():
        continue

    if re.search(r"_Chunk\d+$", directory.name):
        continue

    input_file = directory / "ZZ4lAnalysis.root"

    if not input_file.is_file():
        print(f"WARNING: {input_file} not found!")
        continue

    output_file = outputDir / f"{directory.name}.root"

    print(f"{input_file} -> {output_file}")
    input_file.rename(output_file)
