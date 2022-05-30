import os

def runera(era,lumi):
  os.system("cp samples/2018/ZZTo4la0z"+era+".root samples/2018/ZZTo4l.root")
  os.system("python python/run.py PPZZAnalyzer ZZTo4l -y 2018 -l "+str(lumi))
  os.system("mv results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root results/2018/PPZZAnalyzer_SR4P/ZZTo4la0z"+era+".root")
  os.system("rm samples/2018/ZZTo4l.root")

os.system("make")
runera("D2",10416)
runera("D1",19881)
runera("C",6530)
runera("B2",401)
runera("B1",6384)
runera("A",12104)
os.system("rm results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root")
os.system("hadd results/2018/PPZZAnalyzer_SR4P/ZZTo4l.root results/2018/PPZZAnalyzer_SR4P/ZZTo4la0zA.root results/2018/PPZZAnalyzer_SR4P/ZZTo4la0zB1.root results/2018/PPZZAnalyzer_SR4P/ZZTo4la0zB2.root results/2018/PPZZAnalyzer_SR4P/ZZTo4la0zC.root results/2018/PPZZAnalyzer_SR4P/ZZTo4la0zD1.root results/2018/PPZZAnalyzer_SR4P/ZZTo4la0zD2.root")
