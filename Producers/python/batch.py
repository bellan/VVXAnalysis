#!/bin/env python

# Adapted from /CMGTools/Production/scripts/crisBatch.py

from __future__ import print_function
import sys
import imp
import copy
import os
import shutil
import pickle
import math
from CMGTools.Production.batchmanager import BatchManager
from VVXAnalysis.Producers.datasetToSource import *
#from CMGTools.Production.datasetToSource import *

def chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

def split(comps):
    # import pdb; pdb.set_trace()
    splitComps = []
    for comp in comps:
        chunkSize = 1
        numJobs = 0
### Interpret splitFactor = # of jobs:
#        if hasattr( comp, 'splitFactor') and comp.splitFactor>1:
#            chunkSize = len(comp.files) / comp.splitFactor
#             if len(comp.files) % comp.splitFactor:
#                 chunkSize += 1
### Interpret splitFactor = #files/job
        if hasattr( comp, 'splitFactor') and comp.splitFactor<len(comp.files):
            chunkSize = comp.splitFactor
            for ichunk, chunk in enumerate( chunks( comp.files, chunkSize)):
                numJobs=numJobs+1
                newComp = copy.deepcopy(comp)
                newComp.files = chunk
                newComp.samplename = comp.name
                newComp.name = '{name}_Chunk{index}'.format(name=newComp.name,
                                                       index=ichunk)
                splitComps.append( newComp )
        else:
            numJobs=1
            splitComps.append( comp )
        print(comp.name, ": files=", len(comp.files), ' chunkSize=', chunkSize, "jobs=", numJobs)
    return splitComps


def batchScriptCERN( index, remoteDir=''):
   '''prepare the LSF version of the batch script, to run on LSF'''
#   print "INDEX", index
#   print "remotedir", remoteDir
   script = """#!/bin/tcsh
#BSUB -q 8nh
#BSUB -o job_%J.txt
#ulimit -v 3000000
limit
cat /proc/cpuinfo
cat /proc/meminfo
cd $CMSSW_BASE/src
cmsenv
cd -
echo 'Environment:'
echo
env
echo
echo 'Copying' ${LS_SUBCWD} to ${PWD} 
cp -rf $LS_SUBCWD .
echo '...done'
echo
echo Workdir content:
ls -l
echo
cd `find . -type d | grep /`
pwd
echo 'Running at:' `date`
cmsRun run_cfg.py >& log.txt
set cmsRunStatus=$?
echo 'cmsRun done at: ' `date` with exit status: $cmsRunStatus
if ( $cmsRunStatus != 0 ) echo $cmsRunStatus > exitStatus.txt
gzip log.txt
echo
echo 'ls: '
pwd
ls -l
echo
echo 'Sending the job directory back...'
cp *.root *.txt *.gz $LS_SUBCWD
if ( -z ZZ4lAnalysis.root ) then
 echo 'Empty file:  ZZ4lAnalysis.root'
 if ( -s ZZ4lAnalysis.root ) then
   echo Retrying...
   sleep 10
   cp *.root $LS_SUBCWD
 endif
endif
echo 'destination dir: ls: '
cd $LS_SUBCWD
pwd
ls -l
setenv ROOT_HIST 0
if ( -s ZZ4lAnalysis.root ) then
 root -q -b '${CMSSW_BASE}/src/ZZAnalysis/AnalysisStep/test/prod/rootFileIntegrity.r(\"ZZ4lAnalysis.root\")'
else
 echo moving empty file
 mv ZZ4lAnalysis.root ZZ4lAnalysis.root.empty
endif

echo '...done at' `date`
exit $cmsRunStatus
""" 
   return script


def batchScriptLocal( index ):
   '''prepare a local version of the batch script, to run using nohup'''

   script = """#!/bin/bash
echo 'running'
cmsRun run_cfg.py
""" 
   return script


def tuneProcess(process):
    return;
            
class MyBatchManager( BatchManager ):
   '''Batch manager specific to cmsRun processes.''' 

   def __init__(self):
       self.DefineOptions()

       self.parser_.add_option("-p", "--pdf", dest="PDFstep",
                               help="Step of PDF systematic uncertainty evaluation. It could be 1 or 2.",
                               default=0)
       
       self.parser_.add_option("-i", "--secondary-input-dir", dest="secondaryInputDir",
                               help="Name of the local input directory for your PDF jobs",
                               default=None)



   def setSecondaryInputDir(self):
       if not self.options_.secondaryInputDir == None:
           self.secondaryInputDir_ = os.path.abspath(self.options_.secondaryInputDir) + '/AAAOK/'
       else: 
           self.secondaryInputDir_ = None
       

   def PrepareJob( self, value, dirname=None):
       '''Prepare a job for a given value.
       
       calls PrepareJobUser, which should be overloaded by the user.
       '''
       print('PrepareJob : %s' % value) 
       dname = dirname
       if dname  is None:
           dname = 'Job_{value}'.format( value=value )
       jobDir = '/'.join( [self.outputDir_, dname])
       print('\t',jobDir) 
       self.mkdir( jobDir )
       self.listOfJobs_.append( jobDir )
       if not self.secondaryInputDir_ == None: self.inputPDFDir_ = '/'.join( [self.secondaryInputDir_, dname])

       self.PrepareJobUser( jobDir, value )
       
         
   def PrepareJobUser(self, jobDir, value ):
       '''Prepare one job. This function is called by the base class.'''
#       print value
#       print splitComponents[value]

       # import pdb; pdb.set_trace()
       #prepare the batch script
       scriptFileName = jobDir+'/batchScript.sh'
       scriptFile = open(scriptFileName,'w')
       # storeDir = self.remoteOutputDir_.replace('/castor/cern.ch/cms','')
       mode = self.RunningMode(options.batch)
       if mode == 'LXPLUS':
           scriptFile.write( batchScriptCERN( value ) )
       elif mode == 'LOCAL':
           scriptFile.write( batchScriptLocal( value ) ) 
       scriptFile.close()
       os.system('chmod +x %s' % scriptFileName)

       # Tune can be:
       # for DATA: "DoubleEle", "DoubleMu","MuEG"
       # for MC:   "", "MCAllEvents", "DYJets_B", "DYJets_NoB"
       #
       # setup can be 2011 or 2012 
       
       SAMPLENAME  = splitComponents[value].samplename
       tune  = splitComponents[value].tune
       setup = splitComponents[value].setup
       xsec  = splitComponents[value].xsec
       
       # Set global parameters
       IsMC = True
       PD = ""
       MCFILTER = ""
       SUPERMELA_MASS = 125.6
       if tune != "" :
           tunes = tune.split(",")
           if ("DoubleMu" in tunes) or ("DoubleEle" in tunes) or ("MuEG" in tunes):
               if len(tunes)!=1 : # Do not allow other options for data
                   raise ValueError("invalid tune").with_traceback(tune)
               IsMC = False
               PD = tune
               
           if ("DYJets_B" in tunes) or  ("DYJets_NoB" in tunes) :
               MCFILTER = "HF" # Note: further customization of the filter module is done below           
           if ("MH126" in tunes) :
               SUPERMELA_MASS = 126
           # customization for "MCAllEvents" is done below

           #FIXME: should check tunes for consistency
       XSEC = xsec
       SKIM_REQUIRED = splitComponents[value].doskim
       
       print(SAMPLENAME, "parameters:", tune, IsMC, PD, MCFILTER, SUPERMELA_MASS, XSEC)

       # Read CFG file so that it is customized with the above globals
       namespace = {'IsMC':IsMC, 'PD':PD, 'MCFILTER':MCFILTER, 'SUPERMELA_MASS':SUPERMELA_MASS, 'SAMPLENAME':SAMPLENAME, 'XSEC':XSEC, 'SKIM_REQUIRED':SKIM_REQUIRED}
       exec(compile(open(cfgFileName, "rb").read(), cfgFileName, 'exec'),namespace)
#       handle = open(cfgFileName, 'r')
#       cfo = imp.load_source("pycfg", cfgFileName, handle)
#       process=cfo.process
#       handle.close()                  

       process = namespace.get('process') 
       process.source = splitComponents[value].source
       process.source.fileNames = splitComponents[value].files

       if splitComponents[value].pdfstep < 2:

           # PDF step 1 case: create also a snippet to be used later in step 2 phase
           if splitComponents[value].pdfstep == 1:
               cfgSnippetPDFStep2 = open(jobDir+'/inputForPDFstep2.py','w')
               cfgSnippetPDFStep2.write('process.source.fileNames = ["file:{0:s}/{1:s}"]\n'.format(self.outputDir_+'/AAAOK'+jobDir.replace(self.outputDir_,''), process.weightout.fileName.value()))
               cfgSnippetPDFStep2.write('process.source.secondaryFileNames = [')
               for item in splitComponents[value].files: cfgSnippetPDFStep2.write("'%s',\n" % item)
               cfgSnippetPDFStep2.write(']')
               cfgSnippetPDFStep2.write( '\n' )
               cfgSnippetPDFStep2.close()

       if "MCAllEvents" in tune:
           process.ZZ4muTree.skipEmptyEvents = False
           process.ZZ4eTree.skipEmptyEvents = False
           process.ZZ2e2muTree.skipEmptyEvents = False

       #tuneProcess( process )
       cfgFile = open(jobDir+'/run_cfg.py','w')
#       cfgFile.write('import FWCore.ParameterSet.Config as cms\n\n')
       # cfgFile.write('import os,sys\n')
       cfgFile.write( process.dumpPython() )
       cfgFile.write( '\n' )

       del namespace
       del process

       if splitComponents[value].pdfstep == 2:
           cfgSnippetPDFStep2 = open(self.inputPDFDir_+'/inputForPDFstep2.py','r')
           shutil.copyfileobj(cfgSnippetPDFStep2,cfgFile)
           cfgSnippetPDFStep2.close()
           

       # Tune process
       # JSON for data
       if IsMC==False :
           if setup==2011 :
               cfgFile.write( open('json_2011.py').read() ) 
           elif setup==2012 :
               cfgFile.write( open('json_2012.py').read() )
           else :
                raise ValueError("invalid setup").with_traceback(setup)
        
       if "DYJets_B" in tune :
           cfgFile.write( 'process.HF = cms.Path(process.heavyflavorfilter)\n\n' )
       elif "DYJets_NoB" in tune :
           cfgFile.write( 'process.HF = cms.Path(~process.heavyflavorfilter)\n\n' )
       if "Signal" in tune and not "NoSignal" in tune:
           cfgFile.write( '\nprocess.genCategory0 = process.genCategory =  cms.EDFilter("ZZGenFilterCategory", Topology = cms.int32(SIGNALDEFINITION), ParticleStatus = cms.int32(1), GenJets = cms.InputTag("genJetSel"), src = cms.InputTag("genParticlesPruned"))\n')
           cfgFile.write( 'process.signalFilters += process.genCategory0\n' )
           cfgFile.write( 'process.postSkimSignalCounter = cms.EDProducer("EventCountProducer")\n' )
           cfgFile.write( 'process.signalFilters += process.postSkimSignalCounter\n\n' )
       if "NoSignal" in tune:
           cfgFile.write( '\nprocess.genCategory0 = process.genCategory =  cms.EDFilter("ZZGenFilterCategory", Topology = cms.int32(SIGNALDEFINITION), ParticleStatus = cms.int32(1), GenJets = cms.InputTag("genJetSel"), src = cms.InputTag("genParticlesPruned"))\n')
           cfgFile.write( 'process.signalFilters += ~process.genCategory0\n' )
           cfgFile.write( 'process.postSkimSignalCounter = cms.EDProducer("EventCountProducer")\n' )
           cfgFile.write( 'process.signalFilters += process.postSkimSignalCounter\n\n' )
       if "NoTaus" in tune:
           cfgFile.write( '\nprocess.genTaus = cms.EDFilter("PdgIdAndStatusCandViewSelector", src = cms.InputTag("genParticlesPruned"), pdgId = cms.vint32( 15 ), status = cms.vint32( 3 ))\n')
           cfgFile.write( 'process.genTauCounterFilter =  cms.EDFilter("CandViewCountFilter", src = cms.InputTag("genTaus"), minNumber = cms.uint32(1))\n')
           cfgFile.write( 'process.preSkimCounter  = cms.EDProducer("EventCountProducer")\n')
           cfgFile.write( 'process.signalFilters += process.genTaus\n')
           cfgFile.write( 'process.signalFilters += ~process.genTauCounterFilter\n')
           cfgFile.write( 'process.signalFilters += process.preSkimCounter\n\n')


       cfgFile.close()


class Component(object):

    def __init__(self, name, user, dataset, pattern, splitFactor, tune, xsec, setup, pdfstep, doskim ):
        self.name = name
        print("checking "+self.name)
        self.source = datasetToSource( user, dataset, pattern) # , True for readCache (?)
        self.files = self.source.fileNames # , True for readCache (?)
        self.splitFactor = int(splitFactor)
        self.tune = tune
        self.xsec = float(xsec)
        self.setup = setup
        self.pdfstep = int(pdfstep)
        if self.pdfstep <0 or self.pdfstep>2:
            print("Unknown PDF step", pdfstep)
            sys.exit(1)
        self.doskim = bool(doskim)
        
      
if __name__ == '__main__':
    batchManager = MyBatchManager()

    batchManager.parser_.usage="""
    %prog [options] <cfgFile>

    Run Colin's python analysis system on the batch.
    Job splitting is determined by your configuration file.
    """

    options, args = batchManager.ParseOptions()
    batchManager.setSecondaryInputDir()
    
    cfgFileName = args[0]

    handle = open(cfgFileName, 'r')
    cfo = imp.load_source("pycfg", cfgFileName, handle)
    setup = cfo.LEPTON_SETUP
    #components = [ Component(na, us, da, pattern, sp, tune, xsec, setup) for na, us,da,pattern,sp,tune,xsec in cfo.samples ]

    from ZZAnalysis.AnalysisStep.readSampleInfo import *
    components = []
    sampleDB = readSampleDB('../python/samples_8TeV.csv')
    for sample, settings in sampleDB.items():
        if settings['execute']:
            pdfstep = batchManager.options_.PDFstep
            if pdfstep == 0 or ((not pdfstep == 0) and settings['pdf']):
                components.append(Component(sample, settings['user'], settings['dataset'], settings['pattern'], settings['splitLevel'], settings['tune'],settings['crossSection'], setup, pdfstep, not bool(settings['pdf']))) #settings['pdf'] used here as full sel, without cuts.

    handle.close()


    splitComponents = split( components )

    listOfValues = list(range(0, len(splitComponents)))
    listOfNames = [comp.name for comp in splitComponents]

    batchManager.PrepareJobs( listOfValues, listOfNames )

    waitingTime = 0.05
    batchManager.SubmitJobs( waitingTime )

