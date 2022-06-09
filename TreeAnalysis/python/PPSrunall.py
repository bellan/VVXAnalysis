import os
import sys

def runera(era,lumi):
  os.system("cp samples/2018/ZZTo4la0z"+era+".root samples/2018/ZZTo4l.root")
  os.system("python python/run.py PPZZAnalyzer ZZTo4l -y 2018 -l "+str(lumi))
  os.system("mv results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root results/2018/PPZZAnalyzer_SR4P/ZZTo4la0z"+era+".root")
  os.system("rm samples/2018/ZZTo4l.root")
  
def runfinalstate(state):
  os.system("cp samples/2018/ggTo"+state+".root samples/2018/ZZTo4l.root")
  os.system("python python/run.py PPZZAnalyzer ZZTo4l -y 2018 -l 55700")
  os.system("mv results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root results/2018/PPZZAnalyzer_SR4P/ggTo"+state+".root")
  os.system("rm samples/2018/ZZTo4l.root") 
  
os.system("make")

if(sys.argv[1]=="signal"):
  runera("D2",10416)
  runera("D1",19881)
  runera("C",6530)
  runera("B2",401)
  runera("B1",6384)
  runera("A",12104)
  os.system("hadd results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root results/2018/PPZZAnalyzer_SR4P/ZZTo4la0zA.root results/2018/PPZZAnalyzer_SR4P/ZZTo4la0zB1.root results/2018/PPZZAnalyzer_SR4P/ZZTo4la0zB2.root results/2018/PPZZAnalyzer_SR4P/ZZTo4la0zC.root results/2018/PPZZAnalyzer_SR4P/ZZTo4la0zD1.root results/2018/PPZZAnalyzer_SR4P/ZZTo4la0zD2.root")
  
elif(sys.argv[1]=="ggZZ"):
  runfinalstate("4e")
  runfinalstate("4mu")
  runfinalstate("2e2mu")
  runfinalstate("4tau")
  runfinalstate("2e2tau")
  runfinalstate("2mu2tau")
  os.system("hadd results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root results/2018/PPZZAnalyzer_SR4P/ggTo4e.root results/2018/PPZZAnalyzer_SR4P/ggTo4mu.root results/2018/PPZZAnalyzer_SR4P/ggTo2e2mu.root results/2018/PPZZAnalyzer_SR4P/ggTo4tau.root results/2018/PPZZAnalyzer_SR4P/ggTo2e2tau.root results/2018/PPZZAnalyzer_SR4P/ggTo2mu2tau.root")
