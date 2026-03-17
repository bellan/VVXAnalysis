# How to set up the NanoAOD Analysis
-----------------------------------------------
You need to download three subsystems, with some packages in them [instead of MyProject put the name of you project, for example thesis, or development, ...]

```wget https://raw.githubusercontent.com/cms-sw/cmssw/master/PhysicsTools/NanoAODTools/standalone/checkoutStandalone.sh```
```bash checkoutStandalone.sh -d <MyProject>```
```cd <MyProject>```
```git clone https://github.com/CJLST/ZZAnalysis.git ZZAnalysis```
```(cd ZZAnalysis; git checkout Run3)```
```cd ..```
```git clone https://github.com/bellan/VVXAnalysis.git VVXAnalysis```
```(cd VVXAnalysis; git checkout Run3NanoAOD)```




