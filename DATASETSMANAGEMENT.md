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
- 