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


Recipe for the tree anlysis step on LXPLUS
---------------------------------------------
- prepare the area as described in "Recipe to setup a working LXPLUS area"
- ```cd VVXAnalysis/NanoAnalysis```

The master command is ```./python/run.py <AnalysisName>``` but some options exists:
```
Options:
  -h, --help            show this help message and exit
  -j NJOBS, --jobs=NJOBS
                        Set number of jobs. Default is None, that means it
                        will use all CPU in the systems
  -c, --condor          submit the jobs through Condor
  -e, --eos             use the location of the samples written in the DB
  -s ('YEAR', 'NAME'), --sample=('YEAR', 'NAME')
                        Analyze just this sample (can be repeated and used by
                        Condor jobs
  -f {espresso,microcentury,longlunch,workday,tomorrow,testmatch,nextweek}, --flavour={espresso,microcentury,longlunch,workday,tomorrow,testmatch,nextweek}
                        JobFlavour Condor. Default is (longlunch)
  --dry-run             write Condor files without submitting them
```
to run the analysis, for example VVXAnalyzer, you can prompt

```./python/run.py VVXAnalyzer```

this is the most simple command. It will parse ```configuration/VVXAnalyzer.yaml``` and process the samples therein stored. By default, it will parallelize the job maximizing the number depending on the number of core in the machine. If you want to limit or simply impose a number of core to be used, you can append the ```-j <number>``` option:

```./python/run.py VVXAnalyzer -j 4```

will use 4 core to process all jobs. As soon as a core get available again, the code process a new sample.

The option ```-e``` allows the analysis to run over the sample in eos. The location of the sample is specified in ```data/samples_DB.json```. If you do not use this option, then you must have the samples, or a soft link to them, in a path like ```samples/<year>/<sample name>.root```.

The option ```-s <year> <sample name>``` allows to run over a subset of the samples specified in the yaml configuration. This is particularly useful for testing purposes or for submitting only a subset of samples via Condor (see below). The option is repeatable. The syntax is simple:

```./python/run.py VVXAnalyzer -s 2024 ZZGamma -s 2022 WZGamma```

is the command to run on two samples:  ZZGamma, version 2024, and WZGamma, version 2022. Remind: they must appear in the yaml as well!

The option ```-c``` allows the user to run over [Condor queues](https://batchdocs.web.cern.ch/) instead of on local CPUs. The combined usage with ```-e``` is encouraged.

The option -f <queue flavour> is to be used in combination with -c option only. The flavour of the queue correspond to the maximum duration of the job:
| Name | Duration |
| - | - |
|espresso     | 20 minutes|
|microcentury | 1 hour|
|longlunch    | 2 hours|
|workday      | 8 hours|
|tomorrow     | 1 day|
|testmatch    | 3 days|
|nextweek     | 1 week|


The default is "longlunch".


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
