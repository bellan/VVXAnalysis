import os
import sys

def runera(sample,era,lumi):
  os.system("cp samples/2018/ZZTo4l"+sample+era+".root samples/2018/ZZTo4l.root")
  os.system("python python/run.py PPZZAnalyzer ZZTo4l -y 2018 -l "+str(lumi))
  os.system("mv results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root results/2018/PPZZAnalyzer_SR4P/ZZTo4l"+sample+era+".root")
  os.system("rm samples/2018/ZZTo4l.root")
  
def runera2017(era,lumi):
  os.system("cp samples/2017/"+era+".root samples/2017/ZZTo4l.root")
  os.system("python python/run.py PPZZAnalyzer ZZTo4l -y 2017 -l "+str(lumi))
  os.system("mv results/2017/PPZZAnalyzer_SR4P/ZZTo4l.root results/2017/PPZZAnalyzer_SR4P/ZZTo4l"+era+".root")
  os.system("rm samples/2017/ZZTo4l.root")
  
def runfinalstate(state):
  os.system("cp samples/2018/ggTo"+state+".root samples/2018/ZZTo4l.root")
  os.system("python python/run.py PPZZAnalyzer ZZTo4l -y 2018 -l 55700")
  os.system("mv results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root results/2018/PPZZAnalyzer_SR4P/ggTo"+state+".root")
  os.system("rm samples/2018/ZZTo4l.root") 
  
os.system("make")

if(sys.argv[1]=="a0z" or sys.argv[1]=="aCz" or sys.argv[1]=="mix" or sys.argv[1]=="a0znoPUP"):
  sample=sys.argv[1]
  runera(sample,"D2",10416)
  runera(sample,"D1",19881)
  runera(sample,"C",6530)
  runera(sample,"B2",401)
  runera(sample,"B1",6384)
  runera(sample,"A",12104)
  os.system("hadd results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root results/2018/PPZZAnalyzer_SR4P/ZZTo4l"+sample+"A.root results/2018/PPZZAnalyzer_SR4P/ZZTo4l"+sample+"B1.root results/2018/PPZZAnalyzer_SR4P/ZZTo4l"+sample+"B2.root results/2018/PPZZAnalyzer_SR4P/ZZTo4l"+sample+"C.root results/2018/PPZZAnalyzer_SR4P/ZZTo4l"+sample+"D1.root results/2018/PPZZAnalyzer_SR4P/ZZTo4l"+sample+"D2.root")
  if(sys.argv[1]=="a0znoPUP"): os.system("mv results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root results/2018/PPZZAnalyzer_SR4P/ZZTo4lnoPUP.root")
  else: os.system("cp results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root results/2018/PPZZAnalyzer_SR4P/"+sample+".root")
  
if(sys.argv[1]=="2017"):
  sample=sys.argv[1] 
  runera2017("F3",3820)
  runera2017("F2",8291)
  runera2017("F1",1792)
  runera2017("E",9398)
  runera2017("D",4209)
  runera2017("C2",3407)
  runera2017("C1",5838)
  runera2017("B",2498)
  os.system("hadd results/2017/PPZZAnalyzer_SR4P/ZZTo4l.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lB.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lC1.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lC2.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lD.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lE.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lF1.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lF1.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lF2.root")
  
if(sys.argv[1]=="2017noPUP"):
  sample=sys.argv[1] 
  runera2017("F3noPUP",3820)
  runera2017("F2noPUP",8291)
  runera2017("F1noPUP",1792)
  runera2017("EnoPUP",9398)
  runera2017("DnoPUP",4209)
  runera2017("C2noPUP",3407)
  runera2017("C1noPUP",5838)
  runera2017("BnoPUP",2498)
  os.system("hadd results/2017/PPZZAnalyzer_SR4P/ZZTo4lnoPUP.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lBnoPUP.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lC1noPUP.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lC2noPUP.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lDnoPUP.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lEnoPUP.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lF1noPUP.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lF1noPUP.root results/2017/PPZZAnalyzer_SR4P/ZZTo4lF2noPUP.root")
   
elif(sys.argv[1]=="ggZZ"):
  runfinalstate("4e")
  runfinalstate("4mu")
  runfinalstate("2e2mu")
  runfinalstate("4tau")
  runfinalstate("2e2tau")
  runfinalstate("2mu2tau")
  os.system("hadd results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root results/2018/PPZZAnalyzer_SR4P/ggTo4e.root results/2018/PPZZAnalyzer_SR4P/ggTo4mu.root results/2018/PPZZAnalyzer_SR4P/ggTo2e2mu.root results/2018/PPZZAnalyzer_SR4P/ggTo4tau.root results/2018/PPZZAnalyzer_SR4P/ggTo2e2tau.root results/2018/PPZZAnalyzer_SR4P/ggTo2mu2tau.root")
  
elif(sys.argv[1]=="data"):
  os.system("python python/run.py PPZZAnalyzer EGamma2018A -y 2018")
  os.system("python python/run.py PPZZAnalyzer EGamma2018B -y 2018") 
  os.system("python python/run.py PPZZAnalyzer EGamma2018C -y 2018") 
  os.system("python python/run.py PPZZAnalyzer Egamma2018D -y 2018") 
  os.system("python python/run.py PPZZAnalyzer DoubleMu2018A -y 2018") 
  os.system("python python/run.py PPZZAnalyzer DoubleMu2018B -y 2018") 
  os.system("python python/run.py PPZZAnalyzer DoubleMu2018C -y 2018") 
  os.system("python python/run.py PPZZAnalyzer DoubleMu2018D -y 2018") 
  os.system("python python/run.py PPZZAnalyzer MuonEG2018A -y 2018") 
  os.system("python python/run.py PPZZAnalyzer MuonEG2018B -y 2018") 
  os.system("python python/run.py PPZZAnalyzer MuonEG2018C -y 2018") 
  os.system("python python/run.py PPZZAnalyzer MuonEG2018D -y 2018")
  os.system("hadd results/2018/PPZZAnalyzer_SR4P/DATA.root results/2018/PPZZAnalyzer_SR4P/EGamma2018A.root results/2018/PPZZAnalyzer_SR4P/EGamma2018B.root results/2018/PPZZAnalyzer_SR4P/EGamma2018C.root results/2018/PPZZAnalyzer_SR4P/EGamma2018D.root results/2018/PPZZAnalyzer_SR4P/DoubleMu2018A.root results/2018/PPZZAnalyzer_SR4P/DoubleMu2018B.root results/2018/PPZZAnalyzer_SR4P/DoubleMu2018C.root results/2018/PPZZAnalyzer_SR4P/DoubleMu2018D.root results/2018/PPZZAnalyzer_SR4P/MuonEG2018A.root results/2018/PPZZAnalyzer_SR4P/MuonEG2018B.root results/2018/PPZZAnalyzer_SR4P/MuonEG2018C.root results/2018/PPZZAnalyzer_SR4P/MuonEG2018D.root")
