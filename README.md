Packages for a multi boson final state analysis
-----------------------------------------------
-----------------------------------------------

The basic objects and selections are based on the H --> ZZ --> 4l analysis; this set of packages can be seen as an extension of the H --> ZZ --> 4l analysis code. 
As a matter of fact, the user needs to follow the very same recipe as of the H --> ZZ --> 4l analysis, and on top of that, check-out the code in this repository.

The philosophy is to run the H --> ZZ --> 4l work-flow up to the production of the Z --> ll bosons and deviate from it adding specific objects for the multi boson analyses, e.g., W --> jj object or VVS tag-jets or other vector boson decays. 

The Multi Boson work-flow produces a ROOT tree file, filled with objects like muons, electrons, jets, vector bosons, described by relatively 
simple data formats (step: tree production).
The actual analysis is then performed on the TTrees contained in the ROOR files. For this step I implemented a C++ framework that puts the user in the condition of immediately start an analysis, even on a laptop, and provided the samples are stored locally, to work off-line (step: tree analysis).

The current structured of the repository is:
- ```VVXAnalysis/DataFormats```  --> Data formats for the object written in the TTree, used in both step previously described.
- ```VVXAnalysis/Producers```    --> CMSSW code for TTrees production. This code is not used in the tree analysis step.
- ```VVXAnalysis/TreeAnalysis``` --> Framework for the tree analysis. It is CMSSW independent.
- ```VVXAnalysis/Commons```      --> Library that can be used both in TreeAnalysis and in CMSSW.
- ```VVXAnalysis/MCDataCards```  --> Repository for run cards of specific samples.

The code is located in this repository: https://github.com/bellan/VVXAnalysis.git

Recipe for the tree production step
-----------------------------------

- In a lxplus like environment, setup your area as for H --> ZZ --> 4l analysis, following the recipe in https://github.com/CJLST/ZZAnalysis/tree/Run2_CutBased_UL.
- Check-out the code from this repository.
  - ```git clone https://github.com/bellan/VVXAnalysis.git VVXAnalysis```
- Compile the code with ```scram b -j 8```
- in ```VVXAnalysis/Producers/test/analyzer.py``` there is an example on cmsRun configuration for an interactive run.
- in ```ZZAnalysis/AnalysisStep/test/prod``` there are queue tools useful for submission/check-status/resubmission/merging.
  The main commands are described here:
  - https://github.com/CJLST/ZZAnalysis/blob/master/AnalysisStep/test/prod/PRODUCTION.md 
  - as starting point one can use as template the ```VVXAnalysis/Producers/python/analyzer_VVjj.py``` file.
 - Details about data-sets and their management are reported here: https://github.com/bellan/VVXAnalysis/DATASETSMANAGEMENT.md


Recipe for tree analysis step
-----------------------------

- In an environment with ROOT (with RooFit), CMAKE and boost library installed:
  - check-out the code from this repository as done above.
  - Generate the Makefile.
  - Compile the code.
  - Make a soft link to readSampleInfo python module (hack, to be fixed).
  - Make a soft link to the directory containing the sample (optional).

```
git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis
git clone https://github.com/bellan/VVXAnalysis.git VVXAnalysis
cd VVXAnalysis/TreeAnalysis
cmake CMakeLists.txt
make
cd python/
ln -s ../../../ZZAnalysis/AnalysisStep/python/readSampleInfo.py
```

- In an environment with SCRAM installed:
  - prepare a CMSSW area. 
  - Check-out the code.
  - Compile with scram.
  - Link the bin to bin/ dir.
  - Make a soft link to readSampleInfo python module (hack, to be fixed).
  - Make a soft link to the directory containing the sample (optional).

```
cmsrel CMSSW_X_Y_Z
cd CMSSW_X_Y_Z/src/
cmsenv
git clone https://github.com/bellan/VVXAnalysis.git VVXAnalysis
scram b
cd VVXAnalysis/TreeAnalysis/bin/
ln -s $CMSSW_BASE/bin/$SCRAM_ARCH/eventAnalyzer
cd ../python/
ln -s ../../../ZZAnalysis/AnalysisStep/python/readSampleInfo.py
cd ..
ln -s <samples-location> samples/

```

- To run the code, please use ```./python/run.py``` and follow the instruction therein written. The normal usage is;
  
  ```./python/run.py <name of analysis class> <name of the sample/data set type>```

  As further option, it can take a bool that forces the analysis to grab the cross-section from the CSV file. 
  The default is ```False```, because normally the tree already contains the cross section from the CSV. This option is meant to be used in case of
  a more precise cross section is made available, or a bug in the cross section assignment for a specific sample is found, allowing the user to postpone a new tree production.

Of course, here in this step, it is supposed that you implement something. I need to give you more info, then. As said the code is steered by the ```./python/run.py``` code, that knows 
how to access the samples and their main characteristics, but the actual code is C++ based.
To implement an analysis, you should inherit from the ```EventAnalyzer``` class, that sets up all the relevant branches, the loop over the events and some other useful utilities for
histogramming. The base class has a pure virtual method (```analyze()```) that must be implemented in the concrete class (your analysis). As a matter of fact, all the analysis should be doable
in the ```analyze()``` method (called each event) and in the ```begin()``` and ```end()``` methods, called before and after the loop over the events starts/ends.
Also the ```cut()``` function, called each event, is supposed to be possibly overloaded, as it can hold a pre-selection of the analysis.

Your class needs then to be registered to be ran by the ```eventAnalyzer``` executable. To do that you have to do two things. First, your class must inherit from ```RegistrableAnalysis.h```, so in the inheritance declaration of your class, make sure you have ```RegistrableAnalysis<YourClass>```. Second, in ```AnalysisFactory.cc```, more precisely in the constructor of the class, add a line like 
```Register("YourClass", &RegistrableAnalysis<YourClass>::create);``` 
<strike>To make your code successfully compile on your laptop, you finally need to modify ```CMakeList.txt``` to add the source code of your analysis. Add it to the ```VVXAnalyzer_SRCS``` variable (```scram``` instead does everything by itself).</strike> [No longer needed with the use of GLOB in CMakeLists.txt]
To make more clear the procedure I have put an example (that it is not supposed to be modified) in ```interface/VVXAnalyzer.h``` and ```src/VVXAnalyzer.cc```. 

Note that the histogrammer utility (a member of the ```EventAnalyzer``` class) allows you to fill plots without bothering 
about histograms booking or writing (see some examples in the EventAnalyzer class).


Recipe for unfolding step
-----------------------------

- you can find documentation about RooUnfold here: http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html
```
cd macros/UnfoldingMacros
svn co https://svnsrv.desy.de/public/unfolding/RooUnfold/trunk RooUnfold
cd RooUnfold
cmsenv
make
```

- To test the RooUnfold code try one of the example. 
```
root -l
.x examples/RooUnfoldExample.cxx
```
