Packages for a multi boson final state analysis
-----------------------------------------------
-----------------------------------------------

The basic objects and selections are based on the H->ZZ->4l analysis; this set of packages can be seen as an extension of the H->ZZ->4l analysis code. 
As a matter of fact, the user needs to follow the very same recipe as of the H->ZZ->4l analysis, and on top of that, check-out the code in this repository.

The philosophy is to run the H->ZZ->4l work-flow up to the production of the Z->ll bosons and deviate from it adding specific objects for the multi boson analyses, e.g.,
W->jj object or VVS tag-jets or other vector boson decays. The input of the analysis are patuples produced for H->ZZ->4l and the idea is to collaborate with that team to extend the
pool of samples to fully match our needs too. 

The Multi Boson work-flow then foreseens the production of ROOT tree files, filled with objects like muons, electrons, jets, vector bosons, described by relatively 
simple data formats (step: tree production).
The actual analysis is then performed on TTrees. For this step I implemented a C++ framework that put the user in the condition of immediately start an analysis, even on a laptop, 
and provided the samples are stored locally, to work off-line (step: tree analysis).

The current structured of the repository is:
- ```VVXAnalysis/DataFormats```  --> Data formats for the object written in the TTree, used in both step previously described.
- ```VVXAnalysis/Producers```    --> CMSSW code for TTrees production. This code is not used in the tree analysis step.
- ```VVXAnalysis/TreeAnalysis``` --> Framework for the tree analysis. It is CMSSW independent.

The code is located in this repository: https://github.com/bellan/VVXAnalysis.git

Recipe for the tree production step
-----------------------------------

- In a lxplus like environment, setup your area has for H->ZZ->4l analysis, following the recipe in https://github.com/CJLST/ZZAnalysis.
- I suggest to run ```git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis``` into a ```tmp/``` area and copy the file ```checkout_539.csh```
  into ```$CMSSW_BASE/src``` and from there run ./checkout_539.csh, making sure that no ```ZZAnalysis``` subsytem already exists.
- Check-out the code from this repository.
  - ```git clone https://github.com/bellan/VVXAnalysis.git VVXAnalysis```
- Compile the code with ```scram b```
- in ```VVXAnalysis/Producers/test/analysis_ZZW.py``` there is an example on cmsRun configuration for an interactive run.
- in ```ZZAnalysis/AnalysisStep/test/prod``` there are queue tools useful for submission/check-status/resubmission/merging.
  The main commands are described here:
  - https://github.com/CJLST/ZZAnalysis/blob/master/AnalysisStep/test/prod/PRODUCTION.md 
  - as starting point one can use as template the ```VVXAnalysis/Producers/test/analyzer_ZZW.py``` file.
  - For our analysis, I have modified the batch.py from ZZAnalisys and put our version in ```VVXAnalysis/Producers/python/batch.py```. 
    The sintax and the accepted options are the same, but it uses a CSV file (https://github.com/bellan/VVXAnalysis/Producers/python/samples_8TeV.csv) 
    as Data Base for the samples given in input.
 - Details about datasets and their management are reported here: https://github.com/bellan/VVXAnalysis/DATASETSMANAGEMENT.md


Recipe for tree analysis step
-----------------------------

- In an environment with ROOT and CMAKE installed:
  - check-out the code from this repository as done above.
  - Generate the Makefile.
  - Compile the code.
  - Make a soft link to readSampleInfo python module (hack, to be fixed).
  - Make a soft link to the directory containing the sample (optional).

```
git clone https://github.com/bellan/VVXAnalysis.git VVXAnalysis
cd VVXAnalysis/TreeAnalysis
cmake CMakeLists.txt
make
cd python/
ln -s ../../Producers/python/readSampleInfo.py
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
ln -s $CMSSW_BASE/bin/slc5_amd64_gcc462/<executable file>
cd ../python/
ln -s ../../Producers/python/readSampleInfo.py
cd ..
ln -s <samples-location> samples/

```

- To run the code, please use ```./python/run.py``` and follow the instruction therein written. The normal usage is;
  
  ```./python/run.py <name of executable> <name of the sample/data set type>```

  As further option, it can take a bool that force the analysis to grab the cross-section from the CSV file. 
  The default is ```False```, because normally the tree already contains the cross section from the CSV. This option is meant to be used in case of
  a more precise cross section is made available, or a bug in the cross section assignment for a specific sample is found, allowing the user to postpone a new tree production.

Off course, here in this step, it is supposed that you implement something. I need to give you more info, then. As said the code is steered by the ```./python/run.py``` code, that knows 
how to access the samples and their main characteristics, but the actual code is C++ based.
To implement an analysis, you should inherit from the ```EventAnalyzer``` class, that set up all the relevant branches, the loop over the events and some other useful utilities for
histogramming. The base class has a pure virtual method (```analyze()```) that must be implemented in the concrete class (your analysis). As a matter of fact, all the analysis should be doable
in the ```analyze()``` method (called each event) and in the ```begin()``` and ```end()``` methods, called before and after the loop over the events starts/ends.
Also the ```cut()``` function, called each event, is supposed to be possibly overloaded, as it can holds a pre-selection of the analysis.
Your class needs then to be instantiated into a ```main()``` and the loop function called; namely, this is done into a ```.cpp``` file.
To make your code successfully compiled, you need to modify ```CMakeList.txt``` (to compile on your laptop) and ```bin/BuildFile.xml``` (to compile on lxplus), 
implementing the directive to compile your code, with the proper dependencies.

To make more clear the procedure, I put an example (that it is not supposed to be modified). The example code is in ```bin/vvxAnalysis.cpp``` and ```src/VVXAnalyzer.*```. 
Untill I do not find a better design, you can clone ```bin/vvxAnalysis.cpp``` and modify the copy to instantiate your analyzer instead of ```VVXAnalyzer```.
To compile your code, please follow what done for ```VVXAnalyzer``` in ```CMakeList.txt``` (to compile on your laptop) and in ```bin/BuildFile.xml``` (to compile on lxplus like environment). 
If you are going to run on lxplus, do not forget to link your executable in ```bin/``` as explained above.

Note that the histogrammer utility (a member of the ```EventAnalyzer``` class) allows you to fill plots without bothering 
about histograms booking or writing (see some examples in the EventAnalyzer class).
