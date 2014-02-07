Packages for a multi boson final state analysis
-----------------------------------------------
-----------------------------------------------

The basic objects and selections are based on the H->ZZ->4l analysis; this set of packages can be seen as an extension of the H->ZZ->4l analysis code. 
As a matter of fact, the user needs follow the very same recipe of the H->ZZ->4l analysis setup and on top of that check-out the code in this repository too.

The philosophy is to run the H->ZZ->4l work-flow up to the production of the Z->ll bosons and deviate from it adding specific objects for the multi boson analyses, e.g.,
W->jj object or VVS tag-jets or other vector boson decays. The input of the analysis are patuples produced by the H->ZZ->4l and the idea is to collaborate with them to extend the
pool of samples to fully match our needs too. 

The Multi Boson work-flow then foreseen the production of ROOT tree files, filled with objects like muons, electrons, jets, vector bosons, described by relatively 
simple data formats (step: tree production).
The actual analysis is then performed on the TTrees. For this step I implemented a C++ framework that put the user in the condition of immediately start the analysis that can run on a laptop, 
and once the samples are stored locally, to work off-line (step: tree analysis).

The current structured of the repository is:
VVXAnalysis/DataFormats  --> Data formats for the object written in the TTree, used in both step previously described.
VVXAnalysis/Producers    --> CMSSW code for tree production. This code is not used in the tree analysis step.
VVXAnalysis/TreeAnalysis --> Framework for the tree analysis. It is CMSSW independent.

Recipe for the tree production step
-----------------------------------

- in a lxplus like environment, setup your area has for H->ZZ->4l analysis:
- check-out the code from this repository:
  - git clone https://github.com/bellan/VVXAnalysis.git VVXAnalysis
- Compile the code. It could be that you need to compile using "-k" option in scram, like 
  ```
  scram b -j8 -k 
  ```
  to prevent the compiler to stop with TreeAnalysis code errors (related with include paths not being recognised by scram).
- in VVXAnalysis/Producers/test/analysis_ZZW.py there is an example on cmsRun configuration for an interactive run.
- in ZZAnalysis/AnalysisStep/test/prod there are queue tools useful for submission/check-status/resubmission/merging.
  The main commands are described here:
  - https://github.com/CJLST/ZZAnalysis/blob/master/AnalysisStep/test/prod/PRODUCTION.md 
  - as starting point one can use as template the VVXAnalysis/Producers/test/analyzer_ZZW.py file.
  - ...
- The list of currently patified samples is in:
  https://github.com/CJLST/ZZAnalysis/blob/master/AnalysisStep/test/prod/analyzer_2012.py
  (I am currently planning to convert it into a csv file to keep trace of relevant information).

Recipe for tree analysis step
-----------------------------

- In an environment with ROOT and CMAKE installed, check-out the code from this repository as done above.
- Generate the Makefile.
- Compile the code.

```
git clone https://github.com/bellan/VVXAnalysis.git VVXAnalysis
cd VVXAnalysis/TreeAnalysis
cmake CMakeList.txt
make
```

- To run the code, please use ./python/run.py and follow the instruction therein written.
  ... more info to come ...

Off course, here in this step, it is supposed that you implement something. I need to give you more info, then. The code is steered by the ./python/run.py code, that knows 
about the samples and their main characteristics. The actual code, after the compiling step has been successfully done, is codified in the ./bin/eventAnalyzer executable.
To implement an analysis, you should inherit from the EventAnalyzer class, that set up all the relevant branches, the loop over the events and some useful utilities for
histogramming. The base class has a pure virtual method (analyze()) that must be implemented in the concrete class (your analysis). As a matter of fact, all the analysis should be doable
in the analyze() method (called each event) and in the begin() and end() methods, called before and after the loop over the events starts/ends. Note that the histogrammer utility (a member of the EventAnalyzer class) allows you to fill plots without bothering about histograms booking or writing (see some examples in the EventAnalyzer class).
To make your code successfully compiled, you need to modify the CMakeList.txt file and implement the directive to compile it, with the proper dependencies. Also, you should modify the src/eventAnalyzer.cpp file to instantiate your class, even better if you
implement a new .cpp file with your analysis instance only, in doing that, make sure your new executable is properly compiled by cmake (i.e., you need to modify the 
CMakeList.txt rules) and that the ./python/run.py knows about it.

