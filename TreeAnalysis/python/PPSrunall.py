import os
import sys

def runera(sample,era,lumi):
  os.system("cp samples/2018/ZZTo4l"+sample+era+".root samples/2018/ZZTo4l.root")
  os.system("python python/run.py PPZZAnalyzer ZZTo4l -y 2018 -l "+str(lumi))
  os.system("mv results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root results/2018/PPZZAnalyzer_SR4P/ZZTo4l"+sample+era+".root")
  os.system("rm samples/2018/ZZTo4l.root")
  
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
  
elif(sys.argv[1]=="ggZZ"):
  runfinalstate("4e")
  runfinalstate("4mu")
  runfinalstate("2e2mu")
  runfinalstate("4tau")
  runfinalstate("2e2tau")
  runfinalstate("2mu2tau")
  os.system("hadd results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root results/2018/PPZZAnalyzer_SR4P/ggTo4e.root results/2018/PPZZAnalyzer_SR4P/ggTo4mu.root results/2018/PPZZAnalyzer_SR4P/ggTo2e2mu.root results/2018/PPZZAnalyzer_SR4P/ggTo4tau.root results/2018/PPZZAnalyzer_SR4P/ggTo2e2tau.root results/2018/PPZZAnalyzer_SR4P/ggTo2mu2tau.root")
