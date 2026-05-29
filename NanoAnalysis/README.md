# How to set up the NanoAOD Analysis
-----------------------------------------------
You need to download three subsystems. There is a script that does it for you. Follow the steps here below, instead of <MyProject> put the name of your project, for example thesis, or development, ..., it wiil hust be the work directory 

```
wget https://raw.githubusercontent.com/bellan/VVXAnalysis/refs/heads/Run3NanoAOD/NanoAnalysis/standalone/standalone.sh
bash standalone.sh <MyProject>
```
If you exit the shell, or you start a fresh one, you need to set the env variable. Go inside you Project, then:
```
source ./VVXAnalysis/NanoAnalysis/standalone/init.sh
```

