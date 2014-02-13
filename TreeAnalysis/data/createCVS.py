#! /usr/bin/env python


##################################
## R. Bellan (UNITO) - Feb 2014 ##
##################################


import sys, os, commands, math, csv

samples = [

########
# DATA #
########

    ### DoubleMu
    # ('DoubleMuA', 'cmgtools', '/DoubleMu/Run2012A-22Jan2013-v1/AOD/PAT_CMG_V5_15_0',       'cmgTuple.*root', 18, "DoubleMu"),
    # ('DoubleMuB', 'cmgtools', '/DoubleMuParked/Run2012B-22Jan2013-v1/AOD/PAT_CMG_V5_15_0', 'cmgTuple.*root', 18, "DoubleMu"),
    # ('DoubleMuC', 'cmgtools', '/DoubleMuParked/Run2012C-22Jan2013-v1/AOD/PAT_CMG_V5_15_0', 'cmgTuple.*root', 18, "DoubleMu"),
    # ('DoubleMuD', 'cmgtools', '/DoubleMuParked/Run2012D-22Jan2013-v1/AOD/PAT_CMG_V5_15_0', 'cmgTuple.*root', 15, "DoubleMu"),

    # ### DoubleEle
    # ('DoubleEleA', 'cmgtools', '/DoubleElectron/Run2012A-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 18, "DoubleEle"),
    # ('DoubleEleB', 'cmgtools', '/DoubleElectron/Run2012B-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 18, "DoubleEle"),
    # ('DoubleEleC', 'cmgtools', '/DoubleElectron/Run2012C-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 18, "DoubleEle"),
    # ('DoubleEleD', 'cmgtools', '/DoubleElectron/Run2012D-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 15, "DoubleEle"),

    #  ### MuEG
    # ('MuEGA', 'cmgtools', '/MuEG/Run2012A-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 18, "MuEG"),
    # ('MuEGB', 'cmgtools', '/MuEG/Run2012B-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 18, "MuEG"),
    # ('MuEGC', 'cmgtools', '/MuEG/Run2012C-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 18, "MuEG"),
    # ('MuEGD', 'cmgtools', '/MuEG/Run2012D-22Jan2013-v1/AOD/PAT_CMG_V5_15_0','cmgTuple.*root', 15, "MuEG"),


##############
# Simulation #
##############

    ### ZZ
    ('ZZTo4mu',         'cmgtools', '/ZZTo4mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                    'cmgTuple.*root', 4,  ""),
    ('ZZTo4e',          'cmgtools', '/ZZTo4e_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                     'cmgTuple.*root', 6,  ""),
    ('ZZTo2mu2tau',     'cmgtools', '/ZZTo2mu2tau_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                'cmgTuple.*root', 20, ""),
    ('ZZTo2e2tau',      'cmgtools', '/ZZTo2e2tau_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                 'cmgTuple.*root', 20, ""),
    ('ZZTo2e2mu',       'cmgtools', '/ZZTo2e2mu_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                  'cmgTuple.*root', 12, ""),
    ('ZZTo4tau',        'cmgtools', '/ZZTo4tau_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                   'cmgTuple.*root', 40, ""),
    ('ZZTo4mu_ext',     'cmgtools', '/ZZTo4mu_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',                'cmgTuple.*root', 4,  ""),
    ('ZZTo4e_ext',      'cmgtools', '/ZZTo4e_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',                 'cmgTuple.*root', 6,  ""),
    ('ZZTo2mu2tau_ext', 'cmgtools', '/ZZTo2mu2tau_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',            'cmgTuple.*root', 20, ""),
    ('ZZTo2e2tau_ext',  'cmgtools', '/ZZTo2e2tau_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',             'cmgTuple.*root', 20, ""),
    ('ZZTo2e2mu_ext',   'cmgtools', '/ZZTo2e2mu_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',              'cmgTuple.*root', 12, ""),
    ('ZZTo4tau_ext',    'cmgtools', '/ZZTo4tau_8TeV_ext-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',               'cmgTuple.*root', 40, ""),
    ('ggZZ4l',          'cmgtools', '/GluGluToZZTo4L_8TeV-gg2zz-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',              'cmgTuple.*root', 3,  ""),
    ('ggZZ2l2l',        'cmgtools', '/GluGluToZZTo2L2L_TuneZ2star_8TeV-gg2zz-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 8,  ""),

    ### VBF samples
    ('VBFH126',         'cmgtools', '/VBF_HToZZTo4L_M-126_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',        'cmgTuple.*root', 15, ""),

    # New split ZH/WH/ttH samples
    ('ttH126',    'cmgtools_group', '/TTbarH_HToZZTo4L_M-126_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',            'cmgTuple.*root', 1,  ""),
    ('WH126',     'cmgtools_group', '/WH_HToZZTo4L_M-126_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',                'cmgTuple.*root', 2,  ""),
    ('ZH126',     'cmgtools_group', '/ZH_HToZZTo4L_M-126_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',                'cmgTuple.*root', 2,  ""),

    ### GluGluH 
    ('powheg15H126',         'cmgtools', '/GluGluToHToZZTo4L_M-126_8TeV-powheg15-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0',         'cmgTuple.*root', 2, ""),
    ('powheg15jhuGenV3H126', 'cmgtools', '/SMHiggsToZZTo4L_M-126_8TeV-powheg15-JHUgenV3-pythia6/Summer12_DR53X-PU_S10_START53_V19-v2/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),
    ('minloH126',            'cmgtools', '/GluGluToHToZZTo4L_M-126_8TeV-minloHJJ-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/PAT_CMG_V5_15_0/', 'cmgTuple.*root', 2, ""),

    # REDUCIBLE BG
           
    ### ZToLJ
    ('DYJetsToLLTuneZ2M50',  'cmgtools', '/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 20, ""),
    ('DYJetsToLLTuneZ2M10',  'cmgtools', '/DYJetsToLL_M-10To50filter_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',          'cmgTuple.*root', 20, ""),

    ### ttbar
    ('TTTo2L2Nu2B',  'cmgtools', '/TTTo2L2Nu2B_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',                   'cmgTuple.*root', 6, ""),

    ### Other-bkgs - for dedicated studies
    ('ZZJetsTo4L',   'cmgtools', '/ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0/',       'cmgTuple.*root', 6, ""),
    
    ('WWWJets',      'cmgtools_group', '/WWWJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',        'cmgTuple.*root', 2, ""),
    ('WWZJets',      'cmgtools_group', '/WWZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('WZZJets',      'cmgtools_group', '/WZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('ZZZJets',      'cmgtools_group', '/ZZZNoGstarJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 2, ""),
    ('TTWJets',      'cmgtools_group', '/TTWJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',        'cmgTuple.*root', 2, ""),
    ('TTZJets',      'cmgtools_group', '/TTZJets_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',     'cmgTuple.*root', 2, ""),
    ('TTWWJets',     'cmgtools_group', '/TTWWJets_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',       'cmgTuple.*root', 2, ""),
    ('WZZ_aMCatNLO', 'cmgtools_group', '/WZZ_8TeV-aMCatNLO-herwig/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/PAT_CMG_V5_15_0',     'cmgTuple.*root', 2, ""),
    ('WZ',           'cmgtools_group', '/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0',      'cmgTuple.*root', 8, ""),
    ('WWJets',       'cmgtools_group', '/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/PAT_CMG_V5_15_0', 'cmgTuple.*root', 8, "")
    ]


#read xsection file
filein  = open('Xsection8TeV_v2.txt','r')

#prepare the csv file
fileoutcsv = open('samples_8TeV.csv','w')
fileoutcsv.write("identifier,crossSection = -99.99,totalEvents = -999,luminosity = -99.99,dataset\n")
csvwriter = csv.writer(fileoutcsv) 

#prepare the py file
fileoutpy = open('samples_8TeV.py','w')
fileoutpy.write('samples = [')


for i in range(0,len(samples)-1):
    sample = samples[i][0]
    filein.seek(0)
    foundsampleinfile = 0
    for line in filein:
        line =  line.strip()
        if line.startswith("#"):  continue
        found =  line.find(sample)
        if(found>0):
            foundsampleinfile =+ 1
            spline =  line.split(" ")
            foundsample = False
            comment = False
            newline = []
            for column in spline:
                if column == sample:
                    foundsample = True
                    for j in range(0,6):
                        newline.append(samples[i][j])
                elif column == '#':
                    comment = True
                elif foundsample and not comment and not column == '' and not column == '\n' and not column == '1' and not column == 'all' and not column == '\t\t':
                    newline.append(float(column))
            if len(newline) == 8:
                newline[6] = round(newline[6] * newline.pop(7),10)
            if not len(newline) == 0:
                fileoutpy.write('{},\n'.format(tuple(newline)))
                # for csv file, simplify the output
                lineforcsv = [newline[0],newline[6],"","",newline[2]]
                csvwriter.writerow(lineforcsv)

    if foundsampleinfile == 0:
        print "{0:s} not found!".format(sample)
        fileoutcsv.write("{0:s},,,,{1:s}\n".format(sample, samples[i][2]))
        newline = samples[i] + (-1,)
        fileoutpy.write('{},\n'.format(newline))
    if foundsampleinfile >1:
        print "More than one instance for {0:s} has been found!".format(sample)

fileoutpy.write(']')                  


