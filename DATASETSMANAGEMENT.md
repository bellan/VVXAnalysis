Data sets management and usage of the data set data base
--------------------------------------------------------

For the VVX Analyses, the interesting samples are collected in a CSV file, that acts as data base (together with a set of python tools). The CSV file is located here:
https://github.com/bellan/VVXAnalysis/Producers/python/samples_8TeV.csv.
The samples are a subset of the H --> ZZ --> 4l analysis one. The full list of samples processed with the H --> ZZ --> 4l framework can be find here:
https://github.com/CJLST/ZZAnalysis/blob/master/AnalysisStep/test/prod/analyzer_2012.py.

There is a tool, ```VVXAnalysis/Producers/python/createCSV.py```, that can create a CSV file out of a list of samples in a format like the one in 
https://github.com/CJLST/ZZAnalysis/blob/master/AnalysisStep/test/prod/analyzer_2012.py; since the CSV file has also the information on the cross sections (at (nxN)L0), the tool
parses also https://github.com/CJLST/ZZAnalysis/blob/master/AnalysisStep/test/Plot/Xsection8TeV_v2.txt to extract the proper value. 
Other cross sections are taken from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV.
On a practical side, ```Xsection8TeV_v2.txt``` is copied in the same dir as ```readCSV.py``` and the missing cross sections taken from the above twiki are added before it gets parsed;
also, a clean python fragment file, ```Producers/python/samples_cfi.py``` is used to tell to ```readCSV.py``` which samples to parse.
```createCSV.py``` returns two files: ```samples_8TeV.csv``` and ```samples_8TeV.py```, the latter has the sample list in the same format as the one in
```analyzer_2012.py```, adding at the end of each entry the value of the corresponding cross section. ```samples_8TeV.py``` is not used in these analyses, 
but was meant to be an easy extention of the H --> ZZ --> 4l framework.


The CSV file is used in both steps, trees production and analysis.

Trees production step:
- The csv file is parsed by a python script that make it a data base. Unless otherwise specified in the execute column, a set of job is prepared for each sample. To remove a sample from the
production, just mark it as False in the execute column.
- The crossSection column is used to set the higher order calculated cross section. It is stored inside the externalCrossSection variable. If it is positive, it is the value used to determine the weght of the sample (the weight is determined at the analysis step).

Trees analysis step:
- the same python code used to parse the csv file in the tree production phase is used here to process the samples. 
- one can pass to the ```run.py``` command either the ```identifier```, to analyze a single sample, or the ```process```, to analyze a family of samples.
- The cross section in the file is normally not used, as in principle, if everything went succesfully, it is already stored inside the trees. However, in case a more precise cross section calculation is provided once the trees have been already produced or for other reasons, it is possible to parse the crossSection column at analysis level and use that values just appending ```True``` at the end of the ```run.py``` command.