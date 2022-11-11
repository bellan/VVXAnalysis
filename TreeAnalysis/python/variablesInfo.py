#!/usr/bin/env python2

VarInfo_zz = {"Mass":["m_{4l} [GeV]","m_{4\ell}",10],"Mjj":["m_{jj} [GeV]","m_{JJ}",20],"Z1Mass":["Z1 Mass","m_{2\ell}",10,],"Z2Mass":["Z2 Mass","m_{2\ell}",10,],"Z1lep0_sip":["Z1 lep 0 Sip","Sip",4],"Z1lep0_iso":["Z1 lep 0 Iso","Iso",4],"Z0lep0_pt":["Z1 lep 0 pT","p_{T}",4],"nJets":["N_{jets} (|#eta^{jet}| < 4.7)","N_{jets} (|#eta^{jet}| < 4.7)",1],"nJets_central":["N_{jets} (|#eta^{jet}| < 4.7)","N_{jets} (|#eta^{jet}| < 4.7)",1],"z":["z1","z1",1],"PtJet1":["p_{T}^{jet1} [GeV]","p_{T}^{jet}",1],"EtaJet1":["#eta^{jet1}","#eta^{jet}",9],"PtJet2":["p_{T}^{jet2} [GeV]","p_{T}^{jet}",1],"EtaJet2":["#eta^{jet2}","#eta^{jet}",10],"Z1pt":["Z1 p_{T}","p_{T}",20],"Z2pt":["Z2 p_{T}","p_{T}",10],"Z1z":["Z1 z","z_{Z_{1}}",7],"Z2z":["Z2 z","z_{Z_{2}}",7],"ptJRatio":["","#Sigma p_{T}/# Sum  ",2],"ptRatio":["","#Sum p_{T}",2],"PtZZ":["p_{T}^{4\\ell}","Sum p_{T}",20],"deltaEtaJJ":["|#eta_{jj}|","|#eta_{jj}|",2],"Dphi":["#Delta #phi_{jj}","#Delta #phi_{jj}",10],"Deta":["|#Delta#eta_{jj}|","#Delta #eta_{jj}",5],"Mjj_Central":["m_{jj}","m_{jj}",20],"Deta_Central":["#Delta #eta_{jj}","#Delta #eta_{jj}",5],"Deta2Jet":["#Delta #eta_{jj}, 2 jet","#Delta #eta_{jj} =2 jet",5],"Deta_noCentral":["#Delta #eta_{jj}, >2 jet","#Delta #eta_{jj} > 2 jet",5],"Deta_1noCentral":["#Delta #eta_{jj}, >2 jet","#Delta #eta_{jj} > 2 jet",5],"PtJet1_noCentral":["#eta Jet","#eta^{jet}",9],"EtaJet1_noCentral":["#eta Jet","#eta^{jet}",10]}

VarInfo_vbs = {"Mass":["m_{4\ell}","m_{4\ell}",40],"Mjj":["m_{jj}","m_{JJ}",20],"Z1Mass":["Z1 Mass","m_{2\ell}",10,],"Z2Mass":["Z2 Mass","m_{2\ell}",10,],"Z1lep0_sip":["Z1 lep 0 Sip","Sip",4],"Z1lep0_iso":["Z1 lep 0 Iso","Iso",4],"Z0lep0_pt":["Z1 lep 0 pT","p_{T}",4],"nJets":["# jets","# jets",1],"nJets_central":["# jets","# jets",1],"z":["z1","z1",1],"PtJet1":["pT Jet","p_{T}^{jet}",10],"EtaJet1":["#eta Jet","#eta^{jet}",10],"PtJet2":["pT Jet","p_{T}^{jet}",10],"EtaJet2":["#eta Jet","#eta^{jet}",10],"Z1pt":["Z1 p_{T}","p_{T}",20],"Z2pt":["Z2 p_{T}","p_{T}",10],"Z1z":["Z1 z","z_{Z_{1}}",7],"Z2z":["Z2 z","z_{Z_{2}}",7],"ptJRatio":["","#Sigma p_{T}/# Sum  ",2],"ptRatio":["","#Sum p_{T}",2],"PtZZ":["p_{T}^{4\\ell}","Sum p_{T}",60],"deltaEtaJJ":["|#eta_{jj}|","|#eta_{jj}|",2],"Dphi":["#Delta #phi_{jj}","#Delta #phi_{jj}",10],"Deta":["#Delta #eta_{jj}","#Delta #eta_{jj}",5],"Mjj_Central":["m_{jj}","m_{jj}",20],"Deta_Central":["#Delta #eta_{jj}","#Delta #eta_{jj}",5]}

def getVarInfo_VVX(region):
    VarInfo_VVX = {
        "AAA_cuts"  : {'title':'Cuts', 'unblind':True, 'logy':False},
        'channel_lep':{'title':'lepton flavour', 'unblind':True}
    }
    return VarInfo_VVX

def getVarInfo_VVGamma(region):
    VarInfo_VVGamma = getVarInfo_VVX(region)
    
    # 4L region
    if region in ['SR4P', 'CR3P1F', 'CR2P2F']:
        rebin_mZZG = 1
        if(region in ['SR4P', 'CR3P1F']):
            rebin_mZZG = 2
        elif(region == 'CR2P2F'):
            rebin_mZZG = 2
        VarInfo_VVGamma.update({
            'ZZ_mass' : {'title':'m_{4\ell}'     },
            'Z0_mass' : {'title':'m_{Z0}'        },
            'Z1_mass' : {'title':'m_{Z1}'        },
            'ZZ_pt'   : {'title':'p_{T}^{Z1}'    },
            'Z0_l0_pt': {'title':'p_{T}^{Z0, l0}'},
            'Z0_l1_pt': {'title':'p_{T}^{Z0, l1}'},
            'Z1_l0_pt': {'title':'p_{T}^{Z1, l0}'},
            'Z1_l1_pt': {'title':'p_{T}^{Z1, l1}'},
            'PhFRClosure_PASS_mZZG': {'title':'m_{ZZ#gamma} [GeV/c^{2}]', 'unblind':False, 'rebin':rebin_mZZG}
        })
        
        for name, title in [('4e','4e'), ('2e2m', '2e2\mu'), ('4m', '4\mu')]:
            VarInfo_VVGamma.update({
                "ZZ_mass_"+name : {'title':"m_{%s}"     %(title), 'rebin':1, 'unblind':True},
                "ZZ_pt_"  +name : {'title':"p_{T}^{%s}" %(title), 'rebin':1, 'unblind':True},
            })
        VarInfo_VVGamma.update({
            'ZZ_mass_noPh'        : {'title':'m_{4\ell}\:,\ no\:\gamma', 'rebin':1, 'unblind':True },
            'ZZ_mass_kinPh'       : {'title':'m_{4\ell}\:,\ \gamma\:kin'                       , 'rebin':1, 'unblind':True },
            # Kinematic selection + pixelSeed + electron veto
            'ZZ_mass_veryLoosePh' : {'title':'m_{4\ell}\:,\ \gamma\:loose'                     , 'rebin':1, 'unblind':True },  # Loose = pass 3 cuts
            'ZZ_mass_failPh'      : {'title':'m_{4\ell}\:,\ \gamma\:loose\,\land\:!tight'      , 'rebin':1, 'unblind':True },
            'ZZ_mass_loosePh'     : {'title':'m_{4\ell}\:,\ \gamma\:tight'                     , 'rebin':1, 'unblind':False},  # Tight = cutBasedIDLoose()
            
            'ZZG_mass_kinPh'      : {'title':'m_{4\ell\gamma}\:,\ \gamma\:kin'                 , 'rebin':1, 'unblind':True },
            'ZZG_mass_veryLoosePh': {'title':'m_{4\ell\gamma}\:,\ \gamma\:loose'               , 'rebin':1, 'unblind':True },
            'ZZG_mass_failPh'     : {'title':'m_{4\ell\gamma}\:,\ \gamma\:loose\,\land\:!tight', 'rebin':1, 'unblind':True },
            'ZZG_mass_loosePh'    : {'title':'m_{4\ell\gamma}\:,\ \gamma\:tight'               , 'rebin':1, 'unblind':False}
        })
    
    # 3L region
    elif region in ['SR3P', 'CR001', 'CR010', 'CR011', 'CR100', 'CR101', 'CR110', 'CR000']:
        VarInfo_VVGamma.update({
            "WZ_cutflow": {'title':'Cuts', 'logy':False},
            'ZW_massT': {'title':'mT_{3\ell\\nu}'   , 'rebin':1, 'unblind':True},
            'ZW_pt'   : {'title':'p_{T}^{3\ell\\nu}', 'rebin':1, 'unblind':True},
            'W_l_pt'  : {'title':'p_{t,l10};GeV/c'},
            'lll_mass': {'title':'m_{lll};GeV/c^{2}'},
            # 'paperSel_ZW_massT' : {'title':'m_{T,3l} [GeV/c^{2}]'},
            # 'paperSel_Z_mass'   : {'title':'m_{Z} [GeV/c^{2}]'   },
            # 'paperSel_W_massT'  : {'title':'m_{T,W} [GeV/c^{2}]' },
            # 'paperSel_ZW_pt'    : {'title':'p_{t,ZW} [GeV/c]'    },
            # 'paperSel_Z_l0_pt'  : {'title':'p_{t,l00} [GeV/c]'   },
            # 'paperSel_Z_l1_pt'  : {'title':'p_{t,l01} [GeV/c]'   },
            # 'paperSel_W_l_pt'   : {'title':'p_{t,l10} [GeV/c]'   },
            # 'paperSel_W_MET_pt' : {'title':'p_{t,MET} [GeV/c]'   },
            # 'paperSel_lll_mass' : {'title':'m_{lll} [GeV/c^{2}]' },
            'PhFRClosure_PASS_mWZG': {'title':'m_{WZ#gamma} [GeV/c^{2}]', 'unblind':False}
            # 'debug3L_l1_FRSF': {'title':'FR(l_{1})' }
            # 'debug3L_l2_FRSF': {'title':'FR(l_{2})' }
            # 'debug3L_l3_FRSF': {'title':'FR(l_{3})' }
            # 'debug3L_ZW_FRSF': {'title':'FR(ZW)'    }
        })
        for name, title in [('3e','3e'), ('2e1m', '2e1\mu'), ('2m1e', '2\mu1e'), ('3m', '3\mu')]:
            VarInfo_VVGamma.update({
                'ZW_massT_'+name : {'title':'m_{%s\\nu}'     %(title), 'rebin':1, 'unblind':True},
                'ZW_pt_'   +name : {'title':'p_{T}^{%s\\nu}' %(title), 'rebin':1, 'unblind':True},
                'W_l_pt_'  +name : {'title':'p_{t,l10};GeV/c'},
                'lll_mass_'+name : {'title':'m_{lll};GeV/c^{2}'}
            })
        VarInfo_VVGamma.update({
            'ZW_massT_noPh'   : {'title':'mT_{3\ell\\nu}\:,\ no\:\gamma', 'rebin':1, 'unblind':True },
        })
        VarInfo_VVGamma.update({
            'ZW_massT_kinPh'       : {'title':'mT_{3\ell\\nu}\:,\ \gamma\:kin'                       , 'rebin':1, 'unblind':True },
            'ZW_massT_veryLoosePh' : {'title':'mT_{3\ell\\nu}\:,\ \gamma\:loose'                     , 'rebin':1, 'unblind':True },
            'ZW_massT_failPh'      : {'title':'mT_{3\ell\\nu}\:,\ \gamma\:loose\,\land\:!tight'      , 'rebin':1, 'unblind':True },
            'ZW_massT_loosePh'     : {'title':'mT_{3\ell\\nu}\:,\ \gamma\:tight'                     , 'rebin':1, 'unblind':False},
            
            'ZWG_massT_kinPh'      : {'title':'mT_{3\ell\\nu\gamma}\:,\ \gamma\:kin'                 , 'rebin':1, 'unblind':True },
            'ZWG_massT_veryLoosePh': {'title':'mT_{3\ell\\nu\gamma}\:,\ \gamma\:loose'               , 'rebin':1, 'unblind':True },
            'ZWG_massT_failPh'     : {'title':'mT_{3\ell\\nu\gamma}\:,\ \gamma\:loose\,\land\:!tight', 'rebin':1, 'unblind':True },
            'ZWG_massT_loosePh'    : {'title':'mT_{3\ell\\nu\gamma}\:,\ \gamma\:tight'               , 'rebin':1, 'unblind':False}
        })
        # for name, title in [('e', 'e'), ('m','\mu')]:
        #     VarInfo_VVGamma.update({
        #         'l3_%s_pt'     %(name): {'title': '3^{rd} %s p_{T}'            %(title)},
        #         'l3_%s_Iso'    %(name): {'title': '3^{rd} %s combRelIsoFSRCorr'%(title)},
        #         'l3_%s_pt_MET' %(name): {'title': '3^{rd} %s p_{T}'            %(title)},
        #         'l3_%s_Iso_MET'%(name): {'title': '3^{rd} %s combRelIsoFSRCorr'%(title)}
        #     })
    
    # 2L region
    elif region in ['SR2P', 'SR2P_1L', 'SR2P_1P', 'CR2P_1F']:
        VarInfo_VVGamma.update({
            'AK4_N'         : {'title':'# AK4'   , 'rebin':1, 'unblind':True},
            'AK4_pt'        : {'title':'p_{T}'   , 'rebin':1, 'unblind':True},
            'AK8_N'         : {'title':'# AK8'   , 'rebin':1, 'unblind':True},
            'AK8_pt'        : {'title':'p_{T}'   , 'rebin':1, 'unblind':True, 'logy':True},
            'Z_mass_2e'     : {'title':'m_{2e}'  , 'rebin':1, 'unblind':True},
            'Z_mass_2m'     : {'title':'m_{2\mu}', 'rebin':1, 'unblind':True},
            'VToJ_mass'     : {'title': 'm_{J}'},
            'VTojj_mass'    : {'title': 'm_{jj}'},
            'VToJFake_mass' : {'title': 'm_{J} fake'},
            'VTojjFake_mass': {'title': 'm_{jj} fake'}
            ,
            'Z_mass_noG'    : {'title':'m_{2\ell}\:,\ no\:\gamma'        , 'rebin':1, 'unblind':True },
            'Z_mass_kinPh'  : {'title':'m_{2\ell}\:,\ kin\:\gamma'       , 'rebin':1, 'unblind':True },
            'Z_mass_veryLooseG':{'title':'m_{2\ell}\:,\ \gamma\:loose'},
            'Z_mass_failPh' : {'title':'m_{2\ell}\:,\ \gamma loose\,\land\:!tight', 'rebin':1, 'unblind':True },
            'Z_mass_looseG' : {'title':'m_{2\ell}\:,\ \gamma\:loose'     , 'rebin':1, 'unblind':False},
            
            'ZG_mass_kinG'  : {'title':'m_{2\ell\gamma}\:,\ \gamma\:kin'},
            'ZG_mass_veryLooseG':{'title':'m_{2\ell\gamma}\:,\ \gamma\:loose'},
            'ZG_mass_failG' : {'title':'m_{2\ell\gamma}\:,\ \gamma loose\,\land\:!tight'},
            'ZG_mass_looseG': {'title':'m_{2\ell\gamma}\:,\ \gamma\:tight', 'unblind':False}
        })
        for Vhad in ['VToJ', 'VToJFake']:
            for classifier in ['PNet', 'deepAK8', 'deepAK8MD']:
                for discriminant in ['TvsQCD', 'WvsQCD', 'ZvsQCD']:
                    VarInfo_VVGamma.update({
                        '%s_%s_%s'%(Vhad, classifier, discriminant): {'title': discriminant+' '+classifier, 'rebin':2, 'unblind':True}
                    })
    
    # Photon stuff
    VarInfo_VVGamma.update({
        'kinPhotons_cutflow'   : {'title':'cut'      , 'unblind':True, 'logy':True},
        'kinPhotons_Nm1'       : {'title':'N-1 cuts' , 'unblind':True},
        'kinPhotons_MVA'       : {'title':'MVA score', 'unblind':True}
        # 'kinPhRes_dR'          : {'title':'#DeltaR'  , 'unblind':False},
        # 'noKinPh_all_genPh_N'  : {'title': '# #gamma_{GEN}'     },
        # 'noKinPh_all_genPh_pt' : {'title': '#gamma_{GEN} p_{T}' },
        # 'noKinPh_all_genPh_eta': {'title': '#gamma_{GEN} #eta'  },
        # 'noKinPh_rec_genPh_pt' : {'title': '#gamma_{GEN} p_{T}' },
        # 'noKinPh_rec_genPh_eta': {'title': '#gamma_{GEN} #eta'  }
        ,
        'kinPh_pt'    : {'title':  'p_{T} \gamma_{kin}'  },
        'veryLoosePh_pt'  : {'title':  'p_{T} \gamma_{loose}'},
        'failPh_pt'   : {'title':  'p_{T} \gamma_{fail}' },
        'loosePh_pt'  : {'title':  'p_{T} \gamma_{tight}'}
        ,
        'kinPh_aeta'  : {'title': '|\eta| \gamma_{kin}'  },
        'veryLoosePh_aeta': {'title': '|\eta| \gamma_{loose}'},
        'failPh_aeta' : {'title': '|\eta| \gamma_{fail}' },
        'loosePh_aeta': {'title': '|\eta| \gamma_{tight}'}
    })
    for e in ['EB', 'EE']:
        VarInfo_VVGamma.update({
            'kinPh_sieie_' +e: {'title':'#sigma_{i#etai#eta}', 'rebin':1, 'unblind':True, 'logy':True},
            'kinPh_HoverE_'+e: {'title':'HoverE'             , 'rebin':2, 'unblind':True, 'logy':True},
            'kinPh_chIso_' +e: {'title':'chIso'              , 'rebin':1, 'unblind':True, 'logy':True, 'logx':True},
            'kinPh_neIso_' +e: {'title':'neIso'              , 'rebin':2, 'unblind':True, 'logy':True},
            'kinPh_phIso_' +e: {'title':'phIso'              , 'rebin':2, 'unblind':True, 'logy':True}
        })
    # for name in ['kin', 'loose', 'medium', 'tight']:
    #     VarInfo_VVGamma.update({
    #         'sigmaiEtaiEta_'+name+'Photons': ['#sigma_{i#etai#eta}', 1, True]
    #     })
    VarInfo_VVGamma.update({
        'kinPh_eScale_N'  : {'title':'Number of #gamma passing selection', 'rebin':1, 'unblind':True},
        'kinPhotons_ID': {'title':'#gamma ID'                         , 'rebin':1, 'unblind':True}
    })
    # for name in ['all', 'kin', 'loose']:
    #     VarInfo_VVGamma.update({
    #         'maxG_minL_DR_'+name: {'title':'max_{#gamma}(min_{l}(#DeltaR(#gamma_{%s}, l))' %(name), 'rebin':1, 'unblind':True},
    #         'minL_DR_'     +name: {'title':'min_{l}(#DeltaR(#gamma_{%s}, l)'               %(name), 'rebin':1, 'unblind':True},
    #     })
    
    return VarInfo_VVGamma


def getVariablesInfo(analyzer, region):
    if   analyzer.startswith("ZZ"     ):
        VarInfo = VarInfo_zz
    elif analyzer.startswith("VVX"    ):
        VarInfo = getVarInfo_VVX(region)
    elif analyzer.startswith("VVGamma"):
        VarInfo = getVarInfo_VVGamma(region)
    else                               :
        VarInfo = getVarInfo_vvx(region)
    return VarInfo


if __name__ == '__main__':
    # test: print the varInfo
    from sys import argv
    from json import dumps
    
    if(len(argv) > 2):
        analyzer = argv[1]
        region   = argv[2]
    else:
        analyzer = 'VVGamma'
        region   = 'SR4P'
    
    print "TEST: analyzer =", analyzer, ", region =", region
    VarInfo = getVariablesInfo(analyzer, region)
    print dumps(VarInfo, indent=2)
