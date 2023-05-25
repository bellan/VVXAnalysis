#!/usr/bin/env python2

VarInfo_zz = {"Mass":["m_{4l} [GeV]","m_{4\ell}",10],"Mjj":["m_{jj} [GeV]","m_{JJ}",20],"Z1Mass":["Z1 Mass","m_{2\ell}",10,],"Z2Mass":["Z2 Mass","m_{2\ell}",10,],"Z1lep0_sip":["Z1 lep 0 Sip","Sip",4],"Z1lep0_iso":["Z1 lep 0 Iso","Iso",4],"Z0lep0_pt":["Z1 lep 0 pT","p_{T}",4],"nJets":["N_{jets} (|#eta^{jet}| < 4.7)","N_{jets} (|#eta^{jet}| < 4.7)",1],"nJets_central":["N_{jets} (|#eta^{jet}| < 4.7)","N_{jets} (|#eta^{jet}| < 4.7)",1],"z":["z1","z1",1],"PtJet1":["p_{T}^{jet1} [GeV]","p_{T}^{jet}",1],"EtaJet1":["#eta^{jet1}","#eta^{jet}",9],"PtJet2":["p_{T}^{jet2} [GeV]","p_{T}^{jet}",1],"EtaJet2":["#eta^{jet2}","#eta^{jet}",10],"Z1pt":["Z1 p_{T}","p_{T}",20],"Z2pt":["Z2 p_{T}","p_{T}",10],"Z1z":["Z1 z","z_{Z_{1}}",7],"Z2z":["Z2 z","z_{Z_{2}}",7],"ptJRatio":["","#Sigma p_{T}/# Sum  ",2],"ptRatio":["","#Sum p_{T}",2],"PtZZ":["p_{T}^{4\\ell}","Sum p_{T}",20],"deltaEtaJJ":["|#eta_{jj}|","|#eta_{jj}|",2],"Dphi":["#Delta #phi_{jj}","#Delta #phi_{jj}",10],"Deta":["|#Delta#eta_{jj}|","#Delta #eta_{jj}",5],"Mjj_Central":["m_{jj}","m_{jj}",20],"Deta_Central":["#Delta #eta_{jj}","#Delta #eta_{jj}",5],"Deta2Jet":["#Delta #eta_{jj}, 2 jet","#Delta #eta_{jj} =2 jet",5],"Deta_noCentral":["#Delta #eta_{jj}, >2 jet","#Delta #eta_{jj} > 2 jet",5],"Deta_1noCentral":["#Delta #eta_{jj}, >2 jet","#Delta #eta_{jj} > 2 jet",5],"PtJet1_noCentral":["#eta Jet","#eta^{jet}",9],"EtaJet1_noCentral":["#eta Jet","#eta^{jet}",10]}

VarInfo_vbs = {"Mass":["m_{4\ell}","m_{4\ell}",40],"Mjj":["m_{jj}","m_{JJ}",20],"Z1Mass":["Z1 Mass","m_{2\ell}",10,],"Z2Mass":["Z2 Mass","m_{2\ell}",10,],"Z1lep0_sip":["Z1 lep 0 Sip","Sip",4],"Z1lep0_iso":["Z1 lep 0 Iso","Iso",4],"Z0lep0_pt":["Z1 lep 0 pT","p_{T}",4],"nJets":["# jets","# jets",1],"nJets_central":["# jets","# jets",1],"z":["z1","z1",1],"PtJet1":["pT Jet","p_{T}^{jet}",10],"EtaJet1":["#eta Jet","#eta^{jet}",10],"PtJet2":["pT Jet","p_{T}^{jet}",10],"EtaJet2":["#eta Jet","#eta^{jet}",10],"Z1pt":["Z1 p_{T}","p_{T}",20],"Z2pt":["Z2 p_{T}","p_{T}",10],"Z1z":["Z1 z","z_{Z_{1}}",7],"Z2z":["Z2 z","z_{Z_{2}}",7],"ptJRatio":["","#Sigma p_{T}/# Sum  ",2],"ptRatio":["","#Sum p_{T}",2],"PtZZ":["p_{T}^{4\\ell}","Sum p_{T}",60],"deltaEtaJJ":["|#eta_{jj}|","|#eta_{jj}|",2],"Dphi":["#Delta #phi_{jj}","#Delta #phi_{jj}",10],"Deta":["#Delta #eta_{jj}","#Delta #eta_{jj}",5],"Mjj_Central":["m_{jj}","m_{jj}",20],"Deta_Central":["#Delta #eta_{jj}","#Delta #eta_{jj}",5]}

def getVarInfo_VVX(region):
    VarInfo_VVX = {
        "AAA_cuts"  : {'title':'Cuts', 'text':True, 'unblind':True, 'logy':False},
        'channel_lep':{'title':'lepton flavour', 'text':True, 'unblind':True}
    }
    return VarInfo_VVX

def getVarInfo_VVGamma(region):
    VarInfo_VVGamma = getVarInfo_VVX(region)
    
    # 4L region
    if region in ['SR4P', 'CR3P1F', 'CR2P2F']:
        rebin_mZZG = 1
        channels = (('4e','4e'), ('2e2m', '2e2\mu'), ('4m', '4\mu'))
        if(region in ['SR4P', 'CR2P2F']):
            rebin_mZZG = 1
        elif(region == 'CR3P1F'):
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
            'Z0_dRll' : {},
            'Z1_dRll' : {},
            'PhFRClosure_VLtoL_pt-aeta_data_PASS_mZZG'  : {'title':'m_{ZZ#gamma} [GeV/c^{2}]', 'unblind':False, 'rebin':rebin_mZZG}, #, 'fake_photons':'PhFRClosure_VLtoL_pt-aeta_data_reweighted_mZZG'},
            'PhFRClosure_VLtoL_pt-aeta_dataZG_PASS_mZZG': {'title':'m_{ZZ#gamma} [GeV/c^{2}]', 'unblind':False, 'rebin':rebin_mZZG}, #, 'fake_photons':'PhFRClosure_VLtoL_pt-aeta_dataZG_reweighted_mZZG'},
            'PhFRClosure_KtoVLexcl_pt-aeta_PASS_mZZG'   : {'title':'m_{ZZ#gamma} [GeV/c^{2}]', 'unblind':True , 'rebin':rebin_mZZG}, #, 'fake_photons':'PhFRClosure_KtoVLexcl_pt-aeta_reweighted_mZZG'   },
            'PhFRClosure_VLtoL_pt-aeta_data_FAIL_mZZG'  : {'title':'m_{ZZ#gamma} [GeV/c^{2}]', 'unblind':True , 'rebin':rebin_mZZG},
            'PhFRClosure_VLtoL_pt-aeta_dataZG_FAIL_mZZG': {'title':'m_{ZZ#gamma} [GeV/c^{2}]', 'unblind':True , 'rebin':rebin_mZZG},
            'PhFRClosure_KtoVLexcl_pt-aeta_FAIL_mZZG'   : {'title':'m_{ZZ#gamma} [GeV/c^{2}]', 'unblind':True , 'rebin':rebin_mZZG}
        })
        
        for name, title in channels:
            VarInfo_VVGamma.update({
                "ZZ_mass_"+name : {'title':"m_{%s}"     %(title), 'rebin':1, 'unblind':True},
                "ZZ_pt_"  +name : {'title':"p_{T}^{%s}" %(title), 'rebin':1, 'unblind':True},
            })
        VarInfo_VVGamma.update({
            'ZZ_mass_noPh'        : {'title':'m_{4\ell}\:,\ no\:\gamma', 'rebin':1, 'unblind':True },
            'ZZ_mass_kinPh'       : {'title':'m_{4\ell}\:,\ \gamma\:kin'                       , 'rebin':1, 'unblind':True },
            # Kinematic selection + pixelSeed + electron veto
            'ZZ_mass_veryLoosePh' : {'title':'m_{4\ell}\:,\ \gamma\:loose'                     , 'rebin':1, 'unblind':False},  # Loose = pass 3 cuts
            'ZZ_mass_failPh'      : {'title':'m_{4\ell}\:,\ \gamma\:loose\,\land\:!tight'      , 'rebin':1, 'unblind':True },
            'ZZ_mass_loosePh'     : {'title':'m_{4\ell}\:,\ \gamma\:tight'                     , 'rebin':1, 'unblind':False},  # Tight = cutBasedIDLoose()
            
            'ZZG_mass_kinPh'      : {'title':'m_{4\ell\gamma}\:,\ \gamma\:kin'                 , 'rebin':1, 'split_prompt_ph':True, 'unblind':True },
            'ZZG_mass_veryLoosePh': {'title':'m_{4\ell\gamma}\:,\ \gamma\:loose'               , 'rebin':1, 'split_prompt_ph':True, 'unblind':True },
            'ZZG_mass_failPh'     : {'title':'m_{4\ell\gamma}\:,\ \gamma\:loose\,\land\:!tight', 'rebin':1, 'split_prompt_ph':True, 'unblind':True },
            'ZZG_mass_loosePh'    : {'title':'m_{4\ell\gamma}\:,\ \gamma\:tight'               , 'rebin':1, 'split_prompt_ph':True, 'unblind':False, 'fake_photons': 'ZZG_mass_reweightPh'}
        })
    
    # 3L region
    elif region in ['SR3P', 'CR001', 'CR010', 'CR011', 'CR100', 'CR101', 'CR110', 'CR000']:
        channels = (('3e','3e'), ('2e1m', '2e1\mu'), ('2m1e', '2\mu1e'), ('3m', '3\mu'))
        VarInfo_VVGamma.update({
            "WZ_cutflow": {'title':'Cuts', 'logy':False},
            'ZW_massT': {'title':'mT_{3\ell\\nu}'   , 'rebin':1, 'unblind':True},
            'ZW_pt'   : {'title':'p_{T}^{3\ell\\nu}', 'rebin':1, 'unblind':True},
            'W_l_pt'  : {'title':'p_{t,l10};GeV/c'},
            'lll_mass': {'title':'m_{lll};GeV/c^{2}'},
            'Z_dRll'  : {},
            # 'paperSel_ZW_massT' : {'title':'m_{T,3l} [GeV/c^{2}]'},
            # 'paperSel_Z_mass'   : {'title':'m_{Z} [GeV/c^{2}]'   },
            # 'paperSel_W_massT'  : {'title':'m_{T,W} [GeV/c^{2}]' },
            # 'paperSel_ZW_pt'    : {'title':'p_{t,ZW} [GeV/c]'    },
            # 'paperSel_Z_l0_pt'  : {'title':'p_{t,l00} [GeV/c]'   },
            # 'paperSel_Z_l1_pt'  : {'title':'p_{t,l01} [GeV/c]'   },
            # 'paperSel_W_l_pt'   : {'title':'p_{t,l10} [GeV/c]'   },
            # 'paperSel_W_MET_pt' : {'title':'p_{t,MET} [GeV/c]'   },
            # 'paperSel_lll_mass' : {'title':'m_{lll} [GeV/c^{2}]' },
            'PhFRClosure_LtoT_pt-aeta_PASS_mWZG'     : {'title':'m_{WZ#gamma} [GeV/c^{2}]', 'unblind':False},
            'PhFRClosure_KtoVL_pt-aeta_PASS_mWZG'    : {'title':'m_{WZ#gamma} [GeV/c^{2}]', 'unblind':False},
            'PhFRClosure_KtoVLexcl_pt-aeta_PASS_mWZG': {'title':'m_{WZ#gamma} [GeV/c^{2}]', 'unblind':True },
            'PhFRClosure_LtoT_pt-aeta_FAIL_mWZG'     : {'title':'m_{WZ#gamma} [GeV/c^{2}]', 'unblind':True },
            'PhFRClosure_KtoVL_pt-aeta_FAIL_mWZG'    : {'title':'m_{WZ#gamma} [GeV/c^{2}]', 'unblind':True },
            'PhFRClosure_KtoVLexcl_pt-aeta_FAIL_mWZG': {'title':'m_{WZ#gamma} [GeV/c^{2}]', 'unblind':True }
            # 'debug3L_l1_FRSF': {'title':'FR(l_{1})' }
            # 'debug3L_l2_FRSF': {'title':'FR(l_{2})' }
            # 'debug3L_l3_FRSF': {'title':'FR(l_{3})' }
            # 'debug3L_ZW_FRSF': {'title':'FR(ZW)'    }
        })
        for name, title in channels:
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
        channels = (('2e','2e'), ('2m', '2\mu'))
        VarInfo_VVGamma.update({
            'Z_mass_2e'     : {'title':'m_{2e}'         , 'logy':True},
            'Z_mass_2m'     : {'title':'m_{2\mu}'       , 'logy':True},
            'VToJ_mass'     : {'title':'m_{J}'          , 'logy':True},
            'VTojj_mass'    : {'title':'m_{jj}'         , 'logy':True}
            ,
            'Z_mass_noPh'       : {'title':'m_{2\ell}\:,\ no\:\gamma'                 , 'logy':True},
            'Z_mass_kinPh'      : {'title':'m_{2\ell}\:,\ kin\:\gamma'                , 'logy':True},
            'Z_mass_veryLoosePh': {'title':'m_{2\ell}\:,\ \gamma\:loose'              , 'logy':True},
            'Z_mass_failPh'     : {'title':'m_{2\ell}\:,\ \gamma loose\,\land\:!tight', 'logy':True},
            'Z_mass_loosePh'    : {'title':'m_{2\ell}\:,\ \gamma\:loose'              , 'logy':True}
            ,
            'ZG_mass_kinPh'      : {'title':'m_{2\ell\gamma}\:,\ \gamma\:kin'                , 'logy':True },
            'ZG_mass_veryLoosePh': {'title':'m_{2\ell\gamma}\:,\ \gamma\:loose'              , 'logy':True },
            'ZG_mass_failPh'     : {'title':'m_{2\ell\gamma}\:,\ \gamma loose\,\land\:!tight', 'logy':True },
            'ZG_mass_loosePh'    : {'title':'m_{2\ell\gamma}\:,\ \gamma\:tight'              , 'logy':True, 'unblind':False}
            ,
            'ZjjG_mass_kinPh'      : {'title':'m_{2\elljj\gamma}\:,\ \gamma\:kin'                , 'logy':True },
            'ZjjG_mass_veryLoosePh': {'title':'m_{2\elljj\gamma}\:,\ \gamma\:loose'              , 'logy':True },
            'ZjjG_mass_failPh'     : {'title':'m_{2\elljj\gamma}\:,\ \gamma loose\,\land\:!tight', 'logy':True },
            'ZjjG_mass_loosePh'    : {'title':'m_{2\elljj\gamma}\:,\ \gamma\:tight'              , 'logy':True, 'unblind':False}
            ,
            'ZJG_mass_kinPh'      : {'title':'m_{2\ellJ\gamma}\:,\ \gamma\:kin'                , 'logy':True },
            'ZJG_mass_veryLoosePh': {'title':'m_{2\ellJ\gamma}\:,\ \gamma\:loose'              , 'logy':True },
            'ZJG_mass_failPh'     : {'title':'m_{2\ellJ\gamma}\:,\ \gamma loose\,\land\:!tight', 'logy':True },
            'ZJG_mass_loosePh'    : {'title':'m_{2\ellJ\gamma}\:,\ \gamma\:tight'              , 'logy':True, 'unblind':False}
        })
        for Vhad in ['VToJ']:
            for classifier in ['PNet', 'deepAK8', 'deepAK8MD']:
                for discriminant in ['TvsQCD', 'WvsQCD', 'ZvsQCD']:
                    VarInfo_VVGamma.update({
                        '%s_%s_%s'%(Vhad, classifier, discriminant): {'title': discriminant+' '+classifier, 'rebin':2, 'logy':True, 'unblind':True}
                    })

    elif(region == 'CRLFR'):
        channels = (('2e+e','2e+e'), ('2e+m', '2e+\mu'), ('2m+e', '2\mu+e'), ('2m+m', '2\mu+\mu'))
        VarInfo_VVGamma.update({
            'ZL_mass' :{'title': 'm_{3l} [GeV/c^{2}]'},
            'Z_mass'  :{'title': 'm_{Z} [GeV/c^{2}]' },
            'Z_l0_pt' :{'title': 'p_{t,lZ0} [GeV/c]' },
            'Z_l1_pt' :{'title': 'p_{t,lZ1} [GeV/c]' },
            'Z_dRll'  :{},
            'L_pt'    :{'title': 'p_{t,l3} [GeV/c]'  }
        })
    
    # Photon stuff
    VarInfo_VVGamma.update({
        'kinPhotons_cuts'      : {'title':'cut'      , 'unblind':True, 'logy':True},
        'kinPhotons_Nm1'       : {'title':'N-1 cuts' , 'unblind':True, 'text':True},
        'kinPhoton_MVA'        : {'title':'MVA score', 'unblind':True , 'logy':True},
        'veryLoosePhoton_MVA'  : {'title':'MVA score', 'unblind':False, 'logy':True},
        'loosePhoton_MVA'      : {'title':'MVA score', 'unblind':False, 'logy':True}
        ,
        'cuts_kinPhIso'        : {'logy':True, 'text':True, 'unblind':True },
        'cuts_veryLoosePhIso'  : {'logy':True, 'text':True, 'unblind':False},
        'cuts_failPhIso'       : {'logy':True, 'text':True, 'unblind':True },
        'cuts_loosePhIso'      : {'logy':True, 'text':True, 'unblind':False}
        # ,
        # 'kinPhRes_dR'          : {'title':'#DeltaR'  , 'unblind':False},
        # 'noKinPh_all_genPh_N'  : {'title': '# #gamma_{GEN}'     },
        # 'noKinPh_all_genPh_pt' : {'title': '#gamma_{GEN} p_{T}' },
        # 'noKinPh_all_genPh_eta': {'title': '#gamma_{GEN} #eta'  },
        # 'noKinPh_rec_genPh_pt' : {'title': '#gamma_{GEN} p_{T}' },
        # 'noKinPh_rec_genPh_eta': {'title': '#gamma_{GEN} #eta'  }
        ,
        'lead_fsrPhotons_pt'    : {},
        'lead_fsrPhotons_eta'   : {},
        'lead_fsrPhotons_dRl'   : {},
        'sublead_fsrPhotons_pt' : {},
        'sublead_fsrPhotons_dRl': {},
        'sublead_fsrPhotons_eta': {}
        ,
        'sublead_kin_pt'      : {'title': 'p_{T} \gamma_{kin}^{sublead}'},
        'sublead_veryLoose_pt': {'title': 'p_{T} sublead \gamma_{loose}'},
        'sublead_fail_pt'     : {'title': 'p_{T} sublead \gamma_{fail}' },
        'sublead_loose_pt'    : {'title': 'p_{T} sublead \gamma_{tight}'}
        ,
        'furthestKinPh'   : {},
        'furthestVLPh'    : {},
        'furthestFailPh'  : {},
        'furthestLoosePh' : {'unblind':False}
    })

    for status in ('kin', 'kinVetoVL', 'kinVetoL', 'veryLoose', 'VLchIso', 'VLsieie', 'fail', 'fail3', 'fail4a', 'fail4b', 'loose', 'fsrMatched'):
        variables = [('pt', 'p_{T}'), ('aeta', '|#eta|'), ('dRl', '#DeltaR(l, #gamma)'), ('MVA', 'MVA'), ('chIso', 'chIso'), ('sieie', '#sigma_{i#etai#eta}')]
        unblind = not (status == 'kin' or status == 'veryLoose' or status == 'loose' or status ==  'VLchIso' or status == 'VLsieie')
        if region in ('CR3P1F', 'CR2P2F', 'SR4P'): rebin = 4
        else:                                      rebin = 1
        for varname, vartitle in variables:
            n = 'lead_{}_{}'.format(status, varname)
            d = {'title': '%s #gamma_{%s}^{leading}' %(vartitle, status),
                 'unblind': unblind,
                 'logy': True,
                 'split_prompt_ph': True,
                 'rebin': rebin }
            if status == 'loose':
                d.update({
                    'fake_photons': 'lead_fail_MVA_reweight_data'
                })
            VarInfo_VVGamma.update({n: d})

    for chName, chTitle in channels:
        VarInfo_VVGamma.update({
            'kinPhoton_MVA_'      +chName : {'title':'MVA score', 'unblind':True , 'logy':True},
            'veryLoosePhoton_MVA_'+chName : {'title':'MVA score', 'unblind':False, 'logy':True},
            'loosePhoton_MVA_'    +chName : {'title':'MVA score', 'unblind':False, 'logy':True}
        })

    for e in ['EB', 'EE']:
        VarInfo_VVGamma.update({
            'kinPh_sieie_' +e: {'title':'#sigma_{i#etai#eta}', 'rebin':1, 'unblind':True, 'logy':True},
            'kinPh_HoverE_'+e: {'title':'HoverE'             , 'rebin':2, 'unblind':True, 'logy':True},
            'kinPh_chIso_' +e: {'title':'chIso'              , 'rebin':1, 'unblind':True, 'logy':True},# 'logx':True},
            'kinPh_neIso_' +e: {'title':'neIso'              , 'rebin':2, 'unblind':True, 'logy':True},
            'kinPh_phIso_' +e: {'title':'phIso'              , 'rebin':2, 'unblind':True, 'logy':True}
        })
    # for name in ['kin', 'loose', 'medium', 'tight']:
    #     VarInfo_VVGamma.update({
    #         'sigmaiEtaiEta_'+name+'Photons': ['#sigma_{i#etai#eta}', 1, True]
    #     })
    VarInfo_VVGamma.update({
        'kinPh_central_N'       : {'title':'Number of #gamma_{kin}'    , 'unblind':True , 'logy':True, 'text':True},
        'veryLoosePh_central_N' : {'title':'Number of #gamma_{loose}'  , 'unblind':False, 'logy':True, 'text':True},
        'loosePh_central_N'     : {'title':'Number of #gamma_{tight}'  , 'unblind':False, 'logy':True, 'text':True},
        'kinPh_eScale_N'  : {'title':'Number of #gamma passing selection', 'rebin':1, 'unblind':True, 'text':True},
        'kinPhotons_ID'   : {'title':'#gamma ID'                         , 'rebin':1, 'unblind':True, 'text':True}
    })
    # for name in ['all', 'kin', 'loose']:
    #     VarInfo_VVGamma.update({
    #         'maxG_minL_DR_'+name: {'title':'max_{#gamma}(min_{l}(#DeltaR(#gamma_{%s}, l))' %(name), 'rebin':1, 'unblind':True},
    #         'minL_DR_'     +name: {'title':'min_{l}(#DeltaR(#gamma_{%s}, l)'               %(name), 'rebin':1, 'unblind':True},
    #     })

    # Jet stuff
    VarInfo_VVGamma.update({
        'AK4_N'         : {'title':'# AK4'   , 'rebin':1, 'unblind':True, 'logy':True, 'text':True},
        'AK4_pt'        : {'title':'p_{T}'   , 'rebin':1, 'unblind':True, 'logy':True},
        'AK8_N'         : {'title':'# AK8'   , 'rebin':1, 'unblind':True, 'logy':True, 'text':True},
        'AK8_pt'        : {'title':'p_{T}'   , 'rebin':1, 'unblind':True, 'logy':True}
    })
    
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
