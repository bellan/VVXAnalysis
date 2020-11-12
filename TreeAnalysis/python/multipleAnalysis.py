import os
import time

def testsample(samplename):
  os.system("cp samples/2016/"+samplename+".root samples/2016/test.root")
  os.system("python/run.py VVXnocutsAnalyzer test -r MC")
  os.system("mv results/2016/VVXnocutsAnalyzer_MC/test.root results/2016/VVXnocutsAnalyzer_MC/"+samplename+".root")
  os.system("rm samples/2016/test.root")
  
    
testsample("WZZ_SM")
testsample("WZZ_cW_LI")
testsample("WZZ_cW_QU")
testsample("ZZZ_SM")
testsample("ZZZ_cW_LI")
testsample("ZZZ_cW_QU")  
