# How to set up the NanoAOD Analysis
-----------------------------------------------
The basic objects and selections are based on the H --> ZZ --> 4l analysis; this set of packages can be seen as an extension of the H --> ZZ --> 4l analysis code. 
As a matter of fact, the user needs to follow the very same recipe as of the H --> ZZ --> 4l analysis, and on top of that, check-out the code in this repository.

The philosophy is to run the H --> ZZ --> 4l work-flow up to the production of the Z --> ll bosons and deviate from it adding specific objects for the multi boson analyses, e.g., W --> jj object or VVS tag-jets or other vector boson decays. 

The Multi Boson work-flow produces a ROOT tree file, filled with objects like muons, electrons, jets, vector bosons, using the NanoAOD dataformat e NanoAODTools to read the samples.

Recipe to setup a working LXPLUS area
---------------------------------------------

- log on a lxplus machine (```ssh -Y lxplus.cern.ch``` will do the work)
- you may want to create a directory to contain your project (and entering in it)
- Setup your area as for H --> ZZ --> 4l analysis, following the [recipe](https://github.com/CJLST/ZZAnalysis/tree/Run3). You should limit yourself to [this block](https://github.com/CJLST/ZZAnalysis/tree/Run3#to-install-a-complete-cmssw-area-including-this-package)
- ```cd $CMSSW_BASE/src``` and check-out the code from this repository.
- ```git clone https://github.com/bellan/VVXAnalysis.git VVXAnalysis```
- ```cd VVXAnalysis```
- ```git checkout Run3-NF```
- ```cd $CMSSW_BASE/src```  
- ```scram b -j 8```

Recipe for the tree production step on LXPLUS [TO BE UPDATED]
---------------------------------------------

- in ```ZZAnalysis/AnalysisStep/test/prod``` there are queue tools useful for submission/check-status/resubmission/merging.
  The main commands are described here:
  - https://github.com/CJLST/ZZAnalysis/blob/master/AnalysisStep/test/prod/PRODUCTION.md 
  - as starting point one can use as template the ```VVXAnalysis/Producers/python/analyzer_VVjj.py``` file.
 - Details about data-sets and their management are reported here: https://github.com/bellan/VVXAnalysis/DATASETSMANAGEMENT.md


Recipe for the tree analysis step on LXPLUS
---------------------------------------------
- prepare the area as described in "Recipe to setup a working LXPLUS area"
- ```cd VVXAnalysis/NanoAnalysis```
- make sure the CMSSW environment is active (```cmsenv```)

All the paths used by `run.py` (`configurations/`, `data/`, `samples/`, `results/`, `logs/`) are relative to the current directory: always launch it from `VVXAnalysis/NanoAnalysis`.

### Quick start

```
./python/run.py VVXAnalyzer                        # all the samples in the yaml, all the cores of the machine
./python/run.py VVXAnalyzer -j 4                   # the same, on 4 cores
./python/run.py VVXAnalyzer -s 2024 ZZGamma        # just one sample (useful for tests)
./python/run.py VVXAnalyzer -c -f workday          # on Condor (the samples are read from EOS)
./python/run.py VVXAnalyzer --merge                # when the Condor jobs are done: merge the chunks (see below)
```

The master command is ```./python/run.py <AnalysisName>```, where `<AnalysisName>` is the name of the file `configurations/<AnalysisName>.yaml` (without the extension). The available options are:

```
Options:
  -h, --help            show this help message and exit
  -j NJOBS, --jobs=NJOBS
                        Set number of jobs. Default is None, that means it
                        will use all CPU in the systems
  -c, --condor          submit the jobs through Condor, -e option is
                        automatically applied.
  -e, --eos             use the location of the samples written in the DB
  -s ('YEAR', 'NAME'), --sample=('YEAR', 'NAME')
                        Analyze just this sample (can be repeated and used by
                        Condor jobs
  --chunk=N             Process only this chunk of a sample split in chunks
                        (can be repeated). Requires exactly one -s YEAR NAME.
                        Used by the Condor jobs and to re-run failed chunks
  --split               only if both <name>.root and the directory <name>/
                        exist: use the chunks
  --no-split            only if both <name>.root and the directory <name>/
                        exist: use the single file
  --merge               hadd the chunk outputs of the selected samples (after
                        the jobs are finished)
  --keep-chunks         debug: do not remove the outputs of the single chunks
                        after the merge
  --update-sumw         recompute the total genEventSumw of the chunked MC
                        samples, store it in the DB and exit
  -f {espresso,microcentury,longlunch,workday,tomorrow,testmatch,nextweek}, --flavour={espresso,microcentury,longlunch,workday,tomorrow,testmatch,nextweek}
                        JobFlavour Condor. Default is (longlunch)
  --dry-run             write Condor files without submitting them
```

### What `run.py` needs: the yaml configuration and the samples DB

**The yaml file** `configurations/<AnalysisName>.yaml` says which samples to process and with which analyzer:

```yaml
samples:
  years: [2022, 2023, 2024]
  origin: [MC]            # MC, data, or both
  mode: name              # name, process or all
  by_name:
    - ZZTo4l
    - ZZGamma
    #- WZto3LNu            # commented out: not processed
  by_process:
    - VBS

analysis:
  analyzer: VVXAnalyzer
  regions : [L4P, L3P]
```

- `years` and `origin` always act as filters on the samples DB.
- `mode` decides which list is used: `name` uses `by_name`, `process` uses `by_process`, `all` takes every sample of the DB passing the `years`/`origin` filters. The list that is not selected by `mode` is ignored. A name or a process that does not exist in the DB is an error.
- The samples to be processed are the combinations *sample name x year* found in the DB: if `ZZGamma` exists in 2022, 2023 and 2024 and the three years are listed, three samples are processed.
- `analyzer` is the analyzer class to run (the `analysis_name` of the class in `python/analyses/`).
- `regions` are names of flags defined in `python/Regions.py`. **One output file is produced for each region** (see "Output").

**The samples DB** `data/samples_DB.json` has one block per origin and year, each with the list of its samples:

```json
{"origin": "MC", "year": 2024,
 "samples": [{"process": "ZZ", "name": "ZZTo4l", "eos_path": "/eos/user/...", "version": "", "color": 5, "style": 0}]}
```

`eos_path` is the directory where the sample lives, and it is used with the `-e` option. To add a sample, add an entry in the right block. The entry `genEventSumw` (total sum of weights of the samples split in chunks) is written by the framework: do not edit it by hand, and do not edit the DB while `run.py` is running. NB: when `run.py` writes the DB the file is reformatted (indented with tabs), and a lock file `data/samples_DB.json.lock` appears (it is worth adding it to `.gitignore`).

### Where the input files are

Without `-e` the samples are searched in `samples/<year>/` (you can have there soft links); with `-e` in the `eos_path` of the DB (`-c` implies `-e`). In both cases a sample called `<name>` can be:

| Layout | What `run.py` does |
| - | - |
| a single file `<path>/<name>.root` | one task for the whole sample |
| a directory `<path>/<name>/` containing `<name>_Chunk<N>/ZZ4lAnalysis.root`  | one task for each chunk (see "Samples split in chunks") |
| none of the two | the sample is skipped with a `[SKIPPED]` message |
| both | error, unless you choose with `--split` or `--no-split` |

`--split` (use the chunks) and `--no-split` (use the single file) are needed only in the last case and are ignored when there is no ambiguity. The choice is valid for all the selected samples of that run: if two ambiguous samples need different choices, run them separately with `-s`. With `-c` the choice is automatically passed to the jobs.

### Output

For each sample, region and year:

```
results/<year>/<analyzer>_<region>/<sample name>.root
```

so with `regions: [L4P, L3P]` and the sample `ZZGamma` of 2024 you get `results/2024/VVXAnalyzer_L4P/ZZGamma.root` and `results/2024/VVXAnalyzer_L3P/ZZGamma.root`. The layout is the same whether the sample is a single file or is split in chunks, and whether you ran locally or on Condor. MC histograms are filled with the weight

```
1e3 * luminosity * overallEventWeight / genEventSumw
```

(luminosity in 1/fb, cross sections in pb), data with weight 1.

### Running on the local CPUs

```./python/run.py VVXAnalyzer```

parses the yaml and processes the samples therein stored. By default it parallelizes the job maximizing the number of processes depending on the number of cores in the machine. If you want to limit or simply impose a number of cores to be used, you can append the ```-j <number>``` option:

```./python/run.py VVXAnalyzer -j 4```

will use 4 cores to process all jobs. As soon as a core gets available again, the code processes a new sample (or chunk). Each process needs the memory of one analysis: choose `-j` considering also the RAM, not only the cores. The number of processes is never larger than the number of tasks.

**Logs.** With more than one process, the terminal shows only one line per finished task (`[done]` or `[FAILED]`) and everything else (Python and ROOT output, and the traceback in case of errors) goes in `logs/<analyzer>/<year>/<sample name>.log` (`logs/<analyzer>/<year>/<sample name>/Chunk<N>.log` for chunks). You can follow a task with `tail -f`. With only one process (`-j 1`, or just one task) nothing is redirected: the output goes to the terminal.

The option ```-e``` allows the analysis to run over the sample in eos. The location of the sample is specified in ```data/samples_DB.json```. If you do not use this option, then you must have the samples, or a soft link to them, in a path like ```samples/<year>/<sample name>.root``` (or ```samples/<year>/<sample name>/``` for a sample split in chunks).

The option ```-s <year> <sample name>``` allows to run over a subset of the samples specified in the yaml configuration. This is particularly useful for testing purposes or for submitting only a subset of samples via Condor (see below). The option is repeatable, and both the year and the name must match. The syntax is simple:

```./python/run.py VVXAnalyzer -s 2024 ZZGamma -s 2022 WZGamma```

is the command to run on two samples:  ZZGamma, version 2024, and WZGamma, version 2022. Remind: they must appear in the yaml as well! If none of them is in the yaml selection the program stops with an error, if only some are missing it warns and processes the others.

### Samples split in chunks

Very large samples (tens of GB) are not convenient to be processed as a single file. If the sample is available as a set of chunks, as they come from the production, `run.py` processes each chunk as an independent task (a local process or a Condor job) and merges the results at the end. The layout expected on disk is:

```
<path>/<name>/<name>_Chunk<N>/ZZ4lAnalysis.root
```

**Normalization.** The weight of each MC event is divided by `genEventSumw`, that must be the sum of the weights of *all* the chunks, not of the single chunk, otherwise the merged histograms would be wrong. For this reason the total `genEventSumw` of a sample split in chunks is computed once, reading only the `Runs` tree of every chunk, and stored in `data/samples_DB.json`, in the block of the sample:

```json
"genEventSumw": {"value": 329198063.93507, "nChunks": 137}
```

- It is computed automatically by the first run (local or `-c`) of the sample, when it is not in the DB. It can take a while if there are hundreds of chunks. The jobs never compute it: they just read it.
- It is considered valid as long as the number of chunks on disk is the same. If chunks are added or removed it is recomputed at the next run. If you *replace* chunks keeping the same number, the program cannot notice: use `--update-sumw`.
- ```./python/run.py VVXAnalyzer --update-sumw``` recomputes it for the selected chunked MC samples, stores it and exits.
- It is not needed for data (weight 1).

**Temporary outputs and merge.** Each chunk writes its output in a temporary area, with the same structure of the final output:

```
results/<year>/_chunks/<name>/Chunk<N>/<analyzer>_<region>/<name>.root
```

together with `results/<year>/_chunks/<name>/manifest.json`, that records the analyzer, the regions and the chunks to be expected. The merge:

- checks that the output of every chunk and every region is present, and refuses to merge an incomplete set (it prints the missing chunks);
- merges with `hadd` region by region, in a temporary file renamed only if `hadd` succeeds, producing the usual `results/<year>/<analyzer>_<region>/<name>.root`;
- when all the regions are done, removes `results/<year>/_chunks/<name>`. With `--keep-chunks` the single chunk outputs are kept (debug).

When you run locally the merge is done automatically at the end, for the samples whose chunks were all processed without errors. With Condor you have to ask for it, once the jobs are finished:

```./python/run.py VVXAnalyzer --merge```

It works on the samples selected by the yaml (and by `-s`) that have chunk outputs. Starting again a full run of a chunked sample (without `--chunk`) removes the old temporary outputs of that sample.

**Re-running some chunks.** `--chunk N` (repeatable, requires exactly one `-s year name`) processes only the given chunks, both locally and with `-c`, and never merges by itself:

```
./python/run.py VVXAnalyzer --merge                                   # reports: outputs missing for chunk(s) [5, 17]
./python/run.py VVXAnalyzer -c -s 2024 ZZTo4l --chunk 5 --chunk 17
./python/run.py VVXAnalyzer --merge                                   # when they are done
```

With `--chunk` the total `genEventSumw` is not computed: if it is missing in the DB or stale, the program stops and asks to use `--update-sumw`.

### Running on Condor

The option ```-c``` allows the user to run over [Condor queues](https://batchdocs.web.cern.ch/) instead of on local CPUs. **With ```-c``` the option ```-e``` is applied automatically** (a message in the output reminds it): the jobs read the samples from EOS, and the submission looks at the same files (to find the chunks, to compute `genEventSumw`). Therefore the samples must be on EOS, at the `eos_path` written in the DB.

```./python/run.py VVXAnalyzer -c -f workday```

Run it from the lxplus area where you did `cmsenv`: the jobs recreate the same CMSSW environment. `-j` is ignored. What happens:

1. The samples to process are selected as usual (yaml, `-s`), the total `genEventSumw` of the chunked MC samples is computed if needed, and the files below are written in the current directory.
2. **One Condor job is submitted for each sample, or for each chunk of a sample split in chunks.** The job executes on the worker node the very same `run.py`, restricted to one task (`-s <year> <name> [--chunk N]`).

Generated files (overwritten at every submission):

| File | Content |
| - | - |
| `condor_wrapper.sh` | sets the CMSSW environment and runs `run.py` for one task |
| `samples.txt` | one line per job: year, sample name, chunk (`all` for a single file), log path |
| `analysis.sub` | the Condor submit file |

`--dry-run` writes these files without submitting anything, and without computing `genEventSumw`: use it to check what will be done. The generated files are meant to be rewritten at every submission; as a precaution, do not submit again from the same directory while the jobs of the previous submission are still waiting to start.

The option ```-f <queue flavour>``` (to be used with ```-c``` only) selects the maximum duration of the job:

| Name | Duration |
| - | - |
|espresso     | 20 minutes|
|microcentury | 1 hour|
|longlunch    | 2 hours|
|workday      | 8 hours|
|tomorrow     | 1 day|
|testmatch    | 3 days|
|nextweek     | 1 week|

The default is "longlunch". A job exceeding its flavour is killed: choose it on the longest sample (or chunk). Chunks are short jobs, a sample as single file can take much longer.

**Logs and monitoring.** For each job:

```
logs/<AnalysisName>/<year>/<sample name>.out and .err                 (single file)
logs/<AnalysisName>/<year>/<sample name>/Chunk<N>.out and .err        (chunks)
logs/<AnalysisName>/condor.log                                        (Condor log of all the jobs)
```

```
condor_q                        # status of the jobs
condor_q -hold                  # why a job is held
condor_q -better-analyze <id>   # why a job stays idle
tail -f logs/VVXAnalyzer/2024/ZZGamma.out
```

**Important: the jobs read the yaml and the DB when they start, not when you submit.** Do not modify `configurations/<AnalysisName>.yaml` (samples, regions) nor `data/samples_DB.json` until all the jobs have started: a sample removed from the yaml makes its job stop with an error, and different regions would not match what `--merge` expects. The output is written in `results/`, in the directory you submitted from, which must be visible from the worker nodes.

**When the jobs are finished**, check the logs (a job that failed has no output), and for the samples split in chunks run `--merge` (see "Samples split in chunks"). To redo only what failed: `-c -s <year> <name>` for a single-file sample, adding `--chunk N` (repeated) for a chunked one.

### Typical workflows

Quick test of one sample, locally:
```
./python/run.py VVXAnalyzer -s 2024 ZZGamma -j 1
```

Full production on Condor, with a mix of single files and chunked samples:
```
./python/run.py VVXAnalyzer -c --dry-run        # check the files
./python/run.py VVXAnalyzer -c -f workday       # submit
condor_q                                           # wait...
./python/run.py VVXAnalyzer --merge                # merge the chunks
```

### Troubleshooting

| Message / symptom | Cause and solution |
| - | - |
| `both <name>.root and the directory <name>/ exist` | Choose with `--split` or `--no-split`, or remove one of the two. |
| `[SKIPPED] ... neither ... nor the directory ... exist` | The sample is in the yaml but not at the expected location: wrong path, or (in a local run) you forgot `-e`. |
| `genEventSumw ... is missing in the DB or computed on a different set of chunks` | Happens with `--chunk`. Run `--update-sumw` first (or a run without `--chunk`). |
| `merge NOT done, outputs missing for chunk(s) [...]` | Some jobs failed or are still running. Check their logs, re-run them with `--chunk`, then `--merge` again. |
| `None of the sample in the list is present in the YAML selection` | The year/name given with `-s` is not selected by the yaml (year not in `years`, name commented out, ...). |
| Condor job held | `condor_q -hold` gives the reason. Typical: the log directories or the working directory not accessible. |


Recipe for the tree analysis step on your own laptop
---------------------------------------------

You need to download three subsystems: ```ZZAnalysis```, ```PhysicsTools``` and ```VVXAnalysis```. There is a script that does it for you. Follow the steps here below, instead of <MyProject> put the name of your project, for example thesis, or development, ..., it will be your working directory: 

```
wget https://raw.githubusercontent.com/bellan/VVXAnalysis/refs/heads/Run3NanoAOD/NanoAnalysis/standalone/standalone.sh
bash standalone.sh <MyProject>
```
If you exit the shell, or you start a fresh one, you need to set the env variable. Go inside you Project, then:
```
source ./VVXAnalysis/NanoAnalysis/standalone/init.sh
```
The command to run the analysis works as on LXPLUS, with exception that Condor and eos options are not working (unless the user configured her/his setup accordingly). All the other options are meant to work.
