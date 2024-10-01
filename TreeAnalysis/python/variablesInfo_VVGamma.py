#!/usr/bin/env python

########################################################################
# Dictionary of plots for VVGammaAnalyzer, to be used in variablesInfo #
#                                                                      #
# Author: A. Mecca (alberto.mecca@cern.ch)                             #
########################################################################

from variablesInfo_VVX import getVarInfo_VVX

def getVarInfo_VVGamma(region):
    VarInfo_VVGamma = getVarInfo_VVX(region)
    is_SR = region.startswith('SR')

    # 4L region
    if region in ['SR4P', 'CR3P1F', 'CR2P2F']:
        rebin_mZZG = 1
        channels = (('4e','4e'), ('2e2m', '2e2\mu'), ('4m', '4\mu'))
        if(region in ['SR4P', 'CR2P2F']):
            rebin_mZZG = 1
        elif(region == 'CR3P1F'):
            rebin_mZZG = 2
        VarInfo_VVGamma.update({
            'ZZ_mass' : {'title':'m_{4\ell} [GeV]'     },
            'Z0_mass' : {'title':'m_{Z0} [GeV]'        },
            'Z1_mass' : {'title':'m_{Z1} [GeV]'        },
            'ZZ_pt'   : {'title':'p_{T}^{Z1} [GeV]'    },
            'Z0_l0_pt': {'title':'p_{T}^{Z0, l0} [GeV]', 'xmax': 300.},
            'Z0_l1_pt': {'title':'p_{T}^{Z0, l1} [GeV]', 'xmax': 300.},
            'Z1_l0_pt': {'title':'p_{T}^{Z1, l0} [GeV]', 'xmax': 300.},
            'Z1_l1_pt': {'title':'p_{T}^{Z1, l1} [GeV]', 'xmax': 300.},
            'Z0_l0_eta':{'title':'#eta^{Z0, l0}'       , 'rebin':3, 'scale_ymax': 1.75},
            'Z0_l1_eta':{'title':'#eta^{Z0, l1}'       , 'rebin':3, 'scale_ymax': 1.75},
            'Z1_l0_eta':{'title':'#eta^{Z1, l0}'       , 'rebin':3, 'scale_ymax': 1.75},
            'Z1_l1_eta':{'title':'#eta^{Z1, l1}'       , 'rebin':3, 'scale_ymax': 1.75},
            'Z0_dRll' : {},
            'Z1_dRll' : {},
            'PhFRClosure_VLtoL_pt-aeta_data_PASS_mZZG'  : {'title':'m_{ZZ#gamma} [GeV]', 'unblind':False, 'rebin':rebin_mZZG}, #, 'fake_photons':'PhFRClosure_VLtoL_pt-aeta_data_reweighted_mZZG'},
            'PhFRClosure_VLtoL_pt-aeta_dataZG_PASS_mZZG': {'title':'m_{ZZ#gamma} [GeV]', 'unblind':False, 'rebin':rebin_mZZG}, #, 'fake_photons':'PhFRClosure_VLtoL_pt-aeta_dataZG_reweighted_mZZG'},
            'PhFRClosure_KtoVLexcl_pt-aeta_PASS_mZZG'   : {'title':'m_{ZZ#gamma} [GeV]', 'unblind':True , 'rebin':rebin_mZZG}, #, 'fake_photons':'PhFRClosure_KtoVLexcl_pt-aeta_reweighted_mZZG'   },
            'PhFRClosure_VLtoL_pt-aeta_data_FAIL_mZZG'  : {'title':'m_{ZZ#gamma} [GeV]', 'unblind':True , 'rebin':rebin_mZZG},
            'PhFRClosure_VLtoL_pt-aeta_dataZG_FAIL_mZZG': {'title':'m_{ZZ#gamma} [GeV]', 'unblind':True , 'rebin':rebin_mZZG},
            'PhFRClosure_KtoVLexcl_pt-aeta_FAIL_mZZG'   : {'title':'m_{ZZ#gamma} [GeV]', 'unblind':True , 'rebin':rebin_mZZG}
        })

        for phstatus in ('noph', 'kinVetoL', 'loose'):
            VarInfo_VVGamma.update({
                'Z0_mass_'  +phstatus: {'title':'m_{Z0} [GeV]'        , 'split_prompt_ph': True},
                'Z1_mass_'  +phstatus: {'title':'m_{Z1} [GeV]'        , 'split_prompt_ph': True},
                'ZZ_pt_'    +phstatus: {'title':'p_{T}^{Z1} [GeV]'    , 'split_prompt_ph': True},
                'Z0_l0_pt_' +phstatus: {'title':'p_{T}^{Z0, l0} [GeV]', 'split_prompt_ph': True},
                'Z0_l1_pt_' +phstatus: {'title':'p_{T}^{Z0, l1} [GeV]', 'split_prompt_ph': True},
                'Z1_l0_pt_' +phstatus: {'title':'p_{T}^{Z1, l0} [GeV]', 'split_prompt_ph': True},
                'Z1_l1_pt_' +phstatus: {'title':'p_{T}^{Z1, l1} [GeV]', 'split_prompt_ph': True},
                'Z0_l0_eta_'+phstatus: {'title':'#eta^{Z0, l0}'       , 'split_prompt_ph': True},
                'Z0_l1_eta_'+phstatus: {'title':'#eta^{Z0, l1}'       , 'split_prompt_ph': True},
                'Z1_l0_eta_'+phstatus: {'title':'#eta^{Z1, l0}'       , 'split_prompt_ph': True},
                'Z1_l1_eta_'+phstatus: {'title':'#eta^{Z1, l1}'       , 'split_prompt_ph': True},
            })

        for name, title in channels:
            VarInfo_VVGamma.update({
                "ZZ_mass_"+name : {'title':"m_{%s} [GeV]"     %(title), 'rebin':1, 'unblind':True},
                "ZZ_pt_"  +name : {'title':"p_{T}^{%s} [GeV]" %(title), 'rebin':1, 'unblind':True},
                "Z0_l0_pt_"+name: {'title':"%s - p_{T}^{Z0, l0} [GeV]" %(title), 'rebin':1, 'unblind':True},
                "Z0_l1_pt_"+name: {'title':"%s - p_{T}^{Z0, l1} [GeV]" %(title), 'rebin':1, 'unblind':True},
                "Z1_l0_pt_"+name: {'title':"%s - p_{T}^{Z1, l0} [GeV]" %(title), 'rebin':1, 'unblind':True},
                "Z1_l1_pt_"+name: {'title':"%s - p_{T}^{Z1, l1} [GeV]" %(title), 'rebin':1, 'unblind':True},
                "Z0_l0_eta"+name:{'title':'#eta^{Z0, l0}', 'rebin':3},
                "Z0_l1_eta"+name:{'title':'#eta^{Z0, l1}', 'rebin':3},
                "Z1_l0_eta"+name:{'title':'#eta^{Z1, l0}', 'rebin':3},
                "Z1_l1_eta"+name:{'title':'#eta^{Z1, l1}', 'rebin':3},
            })
        VarInfo_VVGamma.update({
            'ZZ_mass_noPh'        : {'title':'m_{4l}, no #gamma [GeV]', 'rebin':1, 'unblind':True },
            'ZZ_mass_kinPh'       : {'title':'m_{4l}, #gamma_{kin} [GeV]'                      , 'rebin':1, 'unblind':False},
            'ZZ_mass_kinVetoL'    : {'title':'m_{4l}, #gamma_{kin and !Loose} [GeV]'           , 'rebin':1, 'unblind':True },
            # Kinematic selection + pixelSeed + electron veto
            'ZZ_mass_veryLoosePh' : {'title':'m_{4l}, #gamma_{VeryLoose} [GeV]'                , 'rebin':1, 'split_prompt_ph':True, 'unblind':False}, # Loose = pass 3 cuts
            'ZZ_mass_failPh'      : {'title':'m_{4l}, #gamma_{VL and !Loose} [GeV]'            , 'rebin':1, 'split_prompt_ph':True, 'unblind':True },
            'ZZ_mass_loosePh'     : {'title':'m_{4l}, #gamma_{Loose} [GeV]'                    , 'rebin':1, 'split_prompt_ph':True, 'unblind':False}, # Tight = cutBasedIDLoose()

            'ZZG_mass_kinPh'      : {'title':'m_{4l#gamma}, #gamma_{kin [GeV]'                 , 'rebin':1, 'split_prompt_ph':region=='SR4P', 'unblind':False},
            'ZZG_mass_kinVetoL'   : {'title':'m_{4l#gamma}, #gamma_{kin and !Loose [GeV]'      , 'rebin':1, 'split_prompt_ph':region=='SR4P', 'unblind':True },
            'ZZG_mass_veryLoosePh': {'title':'m_{4l#gamma}, #gamma_{VeryLoose [GeV]'           , 'rebin':1, 'split_prompt_ph':region=='SR4P', 'unblind':False},
            'ZZG_mass_failPh'     : {'title':'m_{4l#gamma}, #gamma_{VL and !Loose [GeV]'       , 'rebin':1, 'split_prompt_ph':region=='SR4P', 'unblind':True },
            'ZZG_mass_loosePh'    : {'title':'m_{4l#gamma}, #gamma_{Loose [GeV]'               , 'rebin':1, 'split_prompt_ph':region=='SR4P', 'unblind':False, 'fake_photons': 'ZZG_mass_reweightPh'}
            ,
            'SYS_mZZGwp90_central': {'title':'m_{4l#gamma}, #gamma_{wp90} [GeV]' , 'split_prompt_ph':True, 'split_prompt_ph_pattern': 'SYS_mZZGwp90-%s_central', 'unblind':False},
            'SYS_mZZGwp80_central': {'title':'m_{4l#gamma}, #gamma_{wp80} [GeV]' , 'split_prompt_ph':True, 'split_prompt_ph_pattern': 'SYS_mZZGwp80-%s_central', 'unblind':False},
            'SYS_mZZGloose_central':{'title':'m_{4l#gamma}, #gamma_{Loose} [GeV]', 'split_prompt_ph':True, 'split_prompt_ph_pattern': 'SYS_mZZGloose-%s_central', 'fake_photons': 'SYS_mZZGfailReweight_central','unblind':False},
            'SYS_mZllplusZllGloose_central':{'title':'m_{ll#gamma}+m_{ll}, #gamma_{Loose} [GeV]', 'split_prompt_ph':True, 'split_prompt_ph_pattern': 'SYS_mZllplusZllGloose-%s_central', 'fake_photons': 'SYS_mZllplusZllGfailReweight_central', 'unblind':False}
        })

        VarInfo_VVGamma.update({
            'mZZG_compare': {
                'special':True,
                'unblind':True,
                'title':'m_{4l#gamma} #gamma_{Loose} [GeV]',
                'ratio_title': 'data-driven/MC',
                'rebin': 2,
                'ratio_ymax': 5,
                'data': {
                    'plot' :'ZZG_mass_reweightPh',
                    'legend': 'data-driven'
                },
                'stack':{
                    'plot' :'ZZG_mass_loosePh_nonpro'
                }
            }
        })

    # 3L region
    elif region in ['SR3P', 'CR001', 'CR010', 'CR011', 'CR100', 'CR101', 'CR110', 'CR000']:
        channels = (('3e','3e'), ('2e1m', '2e1\mu'), ('2m1e', '2\mu1e'), ('3m', '3\mu'))
        VarInfo_VVGamma.update({
            "WZ_cutflow": {'title':'Cuts', 'logy':False},
            'ZW_massT': {'title':'m_T^{3lv} [GeV]'   , 'rebin':1, 'unblind':True},
            'ZW_pt'   : {'title':'p_{T}^{3lv} [GeV]', 'rebin':1, 'unblind':True},
            'W_l_pt'  : {'title':'p_{T}^{l_{W}} [GeV]'},
            'lll_mass': {'title':'m_{lll} [GeV]'},
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
            'PhFRClosure_LtoT_pt-aeta_PASS_mWZG'     : {'title':'m_{WZ#gamma} [GeV]', 'unblind':False},
            'PhFRClosure_KtoVL_pt-aeta_PASS_mWZG'    : {'title':'m_{WZ#gamma} [GeV]', 'unblind':False},
            'PhFRClosure_KtoVLexcl_pt-aeta_PASS_mWZG': {'title':'m_{WZ#gamma} [GeV]', 'unblind':True },
            'PhFRClosure_LtoT_pt-aeta_FAIL_mWZG'     : {'title':'m_{WZ#gamma} [GeV]', 'unblind':True },
            'PhFRClosure_KtoVL_pt-aeta_FAIL_mWZG'    : {'title':'m_{WZ#gamma} [GeV]', 'unblind':True },
            'PhFRClosure_KtoVLexcl_pt-aeta_FAIL_mWZG': {'title':'m_{WZ#gamma} [GeV]', 'unblind':True }
            # 'debug3L_l1_FRSF': {'title':'FR(l_{1})' }
            # 'debug3L_l2_FRSF': {'title':'FR(l_{2})' }
            # 'debug3L_l3_FRSF': {'title':'FR(l_{3})' }
            # 'debug3L_ZW_FRSF': {'title':'FR(ZW)'    }
        })
        for name, title in channels:
            VarInfo_VVGamma.update({
                'ZW_massT_'+name : {'title':'m_{%s\\nu}'     %(title), 'rebin':1, 'unblind':True},
                'ZW_pt_'   +name : {'title':'p_{T}^{%s\\nu}' %(title), 'rebin':1, 'unblind':True},
                'W_l_pt_'  +name : {'title':'p_{t,l10} [GeV]'},
                'lll_mass_'+name : {'title':'m_{lll} [GeV]'}
            })
        VarInfo_VVGamma.update({
            'ZW_massT_noPh'   : {'title':'mT_{3lv}, no #gamma [GeV]', 'rebin':1, 'unblind':True },
        })
        VarInfo_VVGamma.update({
            'ZW_massT_kinPh'       : {'title':'mT_{3lv}, #gamma_{kin} [GeV]'                , 'rebin':1, 'unblind':False},
            'ZW_massT_kinVetoL'    : {'title':'mT_{3lv}, #gamma_{kin and !Loose} [GeV]'     , 'rebin':1, 'unblind':True },
            'ZW_massT_veryLoosePh' : {'title':'mT_{3lv}, #gamma_{VeryLoose} [GeV]'          , 'rebin':1, 'unblind':False},
            'ZW_massT_failPh'      : {'title':'mT_{3lv}, #gamma_{VL and !Loose} [GeV]'      , 'rebin':1, 'unblind':True },
            'ZW_massT_loosePh'     : {'title':'mT_{3lv}, #gamma_{Loose} [GeV]'              , 'rebin':1, 'unblind':False},

            'ZWG_massT_kinPh'      : {'title':'mT_{3lv#gamma}, #gamma_{kin} [GeV]'          , 'rebin':1, 'unblind':False},
            'ZWG_massT_kinVetoL'   : {'title':'mT_{3lv#gamma}, #gamma_{kin and !Loose} [GeV]', 'rebin':1, 'unblind':True },
            'ZWG_massT_veryLoosePh': {'title':'mT_{3lv#gamma}, #gamma_{VeryLoose} [GeV]'    , 'rebin':1, 'unblind':False},
            'ZWG_massT_failPh'     : {'title':'mT_{3lv#gamma}, #gamma_{VL and !Loose} [GeV]', 'rebin':1, 'unblind':True },
            'ZWG_massT_loosePh'    : {'title':'mT_{3lv#gamma}, #gamma_{Loose} [GeV]'        , 'rebin':1, 'unblind':False, 'fake_photons': 'ZWG_massT_reweightPh'}
            ,
            'SYS_mWZGwp90_central': {'title':'m_{T}^{3l v #gamma}, #gamma wp90 [GeV]' , 'split_prompt_ph':region=='SR3P', 'split_prompt_ph_pattern': 'SYS_mWZGwp90-%s_central', 'unblind':False},
            'SYS_mWZGwp80_central': {'title':'m_{T}^{3l v #gamma}, #gamma wp80 [GeV]' , 'split_prompt_ph':region=='SR3P', 'split_prompt_ph_pattern': 'SYS_mWZGwp80-%s_central', 'unblind':False},
            'SYS_mWZGloose_central':{'title':'m_{T}^{3l v #gamma}, #gamma Loose [GeV]', 'split_prompt_ph':region=='SR3P', 'split_prompt_ph_pattern': 'SYS_mWZGloose-%s_central', 'fake_photons': 'SYS_mWZGfailReweight_central','unblind':False},
        })
        # for name, title in [('e', 'e'), ('m','\mu')]:
        #     VarInfo_VVGamma.update({
        #         'l3_%s_pt'     %(name): {'title': '3^{rd} %s p_{T}'            %(title)},
        #         'l3_%s_Iso'    %(name): {'title': '3^{rd} %s combRelIsoFSRCorr'%(title)},
        #         'l3_%s_pt_MET' %(name): {'title': '3^{rd} %s p_{T}'            %(title)},
        #         'l3_%s_Iso_MET'%(name): {'title': '3^{rd} %s combRelIsoFSRCorr'%(title)}
        #     })

        VarInfo_VVGamma.update({
            'mTZWG_compare': {
                'special':True,
                'unblind':True,
                'title':'m_{T,3l#gamma} #gamma_{Loose} [GeV]',
                'ratio_title': 'data-driven/MC',
                'rebin': 2,
                'ratio_ymax': 5,
                'data': {
                    'plot' :'ZWG_massT_reweightPh',
                    'legend': 'data-driven'
                },
                'stack':{
                    'plot' :'ZWG_massT_loosePh_nonpro'
                }
            }
        })

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
            'Z_mass_kinPh'      : {'title':'m_{2\ell}\:,\ \gamma\:kin'                , 'logy':True, 'unblind':False},
            'Z_mass_kinVetoL'   : {'title':'m_{2\ell}\:,\ \gamma\:kin\,\land\:!tight' , 'logy':True},
            'Z_mass_veryLoosePh': {'title':'m_{2\ell}\:,\ \gamma\:loose'              , 'logy':True, 'unblind':False},
            'Z_mass_failPh'     : {'title':'m_{2\ell}\:,\ \gamma loose\,\land\:!tight', 'logy':True},
            'Z_mass_loosePh'    : {'title':'m_{2\ell}\:,\ \gamma\:loose'              , 'logy':True, 'unblind':False}
            ,
            'ZG_mass_kinPh'      : {'title':'m_{2\ell\gamma}\:,\ \gamma\:kin'                , 'logy':True, 'split_prompt_ph':True, 'unblind':False},
            'ZG_mass_kinVetoL'   : {'title':'m_{2\ell\gamma}\:,\ \gamma\:kin,\land\:!tight'  , 'logy':True, 'split_prompt_ph':True },
            'ZG_mass_veryLoosePh': {'title':'m_{2\ell\gamma}\:,\ \gamma\:loose'              , 'logy':True, 'split_prompt_ph':True, 'unblind':False},
            'ZG_mass_failPh'     : {'title':'m_{2\ell\gamma}\:,\ \gamma loose\,\land\:!tight', 'logy':True, 'split_prompt_ph':True },
            'ZG_mass_loosePh'    : {'title':'m_{2\ell\gamma}\:,\ \gamma\:tight'              , 'logy':True, 'split_prompt_ph':True, 'unblind':False, 'fake_photons': 'ZG_mass_reweightPh'}
            ,
            'ZjjG_mass_kinPh'      : {'title':'m_{2\elljj\gamma}\:,\ \gamma\:kin'                , 'logy':True, 'split_prompt_ph':True, 'unblind':False},
            'ZjjG_mass_kinVetoL'   : {'title':'m_{2\elljj\gamma}\:,\ \gamma\:kin\,\land\:!tight' , 'logy':True, 'split_prompt_ph':True },
            'ZjjG_mass_veryLoosePh': {'title':'m_{2\elljj\gamma}\:,\ \gamma\:loose'              , 'logy':True, 'split_prompt_ph':True, 'unblind':False},
            'ZjjG_mass_failPh'     : {'title':'m_{2\elljj\gamma}\:,\ \gamma loose\,\land\:!tight', 'logy':True, 'split_prompt_ph':True },
            'ZjjG_mass_loosePh'    : {'title':'m_{2\elljj\gamma}\:,\ \gamma\:tight'              , 'logy':True, 'split_prompt_ph':True, 'unblind':False, 'fake_photons': 'ZjjG_mass_reweightPh'}
            ,
            'ZJG_mass_kinPh'      : {'title':'m_{2\ellJ\gamma}\:,\ \gamma\:kin'                , 'logy':True, 'split_prompt_ph':True, 'unblind':False},
            'ZJG_mass_kinVetoL'   : {'title':'m_{2\ellJ\gamma}\:,\ \gamma\:kin\,\land\:!tight' , 'logy':True, 'split_prompt_ph':True },
            'ZJG_mass_veryLoosePh': {'title':'m_{2\ellJ\gamma}\:,\ \gamma\:loose'              , 'logy':True, 'split_prompt_ph':True, 'unblind':False},
            'ZJG_mass_failPh'     : {'title':'m_{2\ellJ\gamma}\:,\ \gamma loose\,\land\:!tight', 'logy':True, 'split_prompt_ph':True },
            'ZJG_mass_loosePh'    : {'title':'m_{2\ellJ\gamma}\:,\ \gamma\:tight'              , 'logy':True, 'split_prompt_ph':True, 'unblind':False, 'fake_photons': 'ZJG_mass_reweightPh'}
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
            'ZL_mass' :{'title': 'm_{3l} [GeV]', 'split_prompt_ph':True},
            'Z_mass'  :{'title': 'm_{Z} [GeV]' , 'split_prompt_ph':True},
            'Z_l0_pt' :{'title': 'p_{T}^{lZ0} [GeV]' , 'split_prompt_ph':True},
            'Z_l1_pt' :{'title': 'p_{T}^{lZ1} [GeV]' , 'split_prompt_ph':True},
            'Z_dRll'  :{'split_prompt_ph':True},
            'L_pt'    :{'title': 'p_{T}^{l3} [GeV]'  , 'split_prompt_ph':True},
            'MET_fine':{'title': 'MET [GeV]'         , 'split_prompt_ph':True},
        })

    # Photon stuff
    VarInfo_VVGamma.update({
        'kinPhotons_cuts'      : {'title':'cut'      , 'unblind':True , 'logy':True},
        'kinPhotons_Nm1'       : {'title':'N-1 cuts' , 'unblind':False, 'text':True},
        'kinPhoton_MVA'        : {'title':'MVA score', 'unblind':False, 'logy':True, 'rebin':2},
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
        'lead_fsrPhotons_pt'    : {'title':'p_{T}^{#gamma} [GeV]', 'rebin':5, 'split_prompt_ph':True},
        'lead_fsrPhotons_aeta'  : {'title':'|#eta^{#gamma}|'     , 'rebin':5, 'split_prompt_ph':True},
        'lead_fsrPhotons_dRl'   : {'title':'#DeltaR(#gamma, l)'  , 'rebin':5, 'split_prompt_ph':True},
        'sublead_fsrPhotons_pt' : {'split_prompt_ph':True},
        'sublead_fsrPhotons_dRl': {'split_prompt_ph':True},
        'sublead_fsrPhotons_aeta':{'split_prompt_ph':True}
        ,
        'sublead_kinVetoL_pt' : {'title': 'p_{T} \gamma_{kin}^{sublead}'},
        'sublead_fail_pt'     : {'title': 'p_{T} sublead \gamma_{fail}' },
        'sublead_loose_pt'    : {'title': 'p_{T} sublead \gamma_{tight}'}
        ,
        'furthestKinPh'   : {},
        'furthestVLPh'    : {},
        'furthestFailPh'  : {},
        'furthestLoosePh' : {'unblind':False}
        ,
        'SYS_MVAcut_central'  : {'title':'MVA cut passed'                  , 'split_prompt_ph':is_SR          , 'split_prompt_ph_pattern': 'SYS_MVAcut-%s_central'  , 'unblind':False, 'logy':True, 'ymin': 1},
        'SYS_wp90pt_central'  : {'title':'p_{T} #gamma_{wp90} [GeV]'       , 'split_prompt_ph':is_SR          , 'split_prompt_ph_pattern': 'SYS_wp90pt-%s_central'  , 'unblind':False},
        'SYS_wp80pt_central'  : {'title':'p_{T} #gamma_{wp80} [GeV]'       , 'split_prompt_ph':is_SR          , 'split_prompt_ph_pattern': 'SYS_wp80pt-%s_central'  , 'unblind':False},
        'SYS_loosept_central' : {'title':'p_{T} #gamma_{Loose} [GeV]'      , 'split_prompt_ph':is_SR          , 'split_prompt_ph_pattern': 'SYS_loosept-%s_central' , 'unblind':False},
    })

    for status in ('kinVetoL', 'fail', 'fail3', 'fail4a', 'fail4b', 'loose', 'fsrMatched', 'FSRkin', 'FSRloose', 'wp90', 'wp80', '90not80'):
        variables = [('pt', 'p_{T}', 'GeV'), ('aeta', '|#eta|', ''), ('dRl', '#DeltaR(l, #gamma)', '')
                     , ('MVA', 'MVA', ''), ('chIso', 'chIso', 'GeV'), ('sieie', '#sigma_{i#etai#eta}', '')
                     , ('pt_fine', 'p_{T}', 'GeV'), ('aeta_fine', '|#eta|', '')
                     ]
        unblind = status not in ('loose', 'wp90', 'wp80')
        for varname, vartitle, udm in variables:
            if region in ('CR3P1F', 'CR2P2F', 'SR4P') and not varname in ('pt', 'aeta', 'pt_fine', 'aeta_fine'):
                rebin = 4
            else:
                rebin = 1

            n = 'lead_{}_{}'.format(status, varname)
            title = '%s #gamma_{%s}^{leading}' %(vartitle, status)
            if(udm is not None and len(udm) > 0):
                title += ' [%s]'%(udm)
            d = {'title': title,
                 'unblind': unblind,
                 'logy': varname in ('MVA',),
                 'split_prompt_ph': True,
                 'rebin': rebin }
            if status == 'loose':
                d.update({
                    'fake_photons': 'lead_fail_{var}_reweight_data'.format(var=varname)
                })
            if('aeta' in varname):
                d.update({'scale_ymax': 1.8})
            if varname in ('dRl','chIso'):
                d.update({'draw_overflow': True})
            VarInfo_VVGamma.update({n: d})

    for chName, chTitle in channels:
        VarInfo_VVGamma.update({
            'kinPhoton_MVA_'      +chName : {'title':'MVA score', 'unblind':True , 'logy':True},
            'veryLoosePhoton_MVA_'+chName : {'title':'MVA score', 'unblind':False, 'logy':True},
            'loosePhoton_MVA_'    +chName : {'title':'MVA score', 'unblind':False, 'logy':True}
        })

    for e in ['EB', 'EE']:
        VarInfo_VVGamma.update({
            'kinPh_sieie_' +e: {'title':'#sigma_{i#etai#eta}', 'rebin':1, 'split_prompt_ph':True, 'unblind':True, 'logy':False},
            'kinPh_HoverE_'+e: {'title':'H/E'                , 'rebin':3, 'split_prompt_ph':False, 'unblind':True, 'logy':False, 'ymin':0},
            'kinPh_chIso_' +e: {'title':'I_{ch} [GeV]'       , 'rebin':3, 'split_prompt_ph':False, 'unblind':True, 'logy':False, 'ymin':0},# 'logx':True},
            'kinPh_neIso_' +e: {'title':'I_{n} [GeV]'        , 'rebin':3, 'split_prompt_ph':False, 'unblind':True, 'logy':False, 'ymin':0},
            'kinPh_phIso_' +e: {'title':'I_{#gamma} [GeV]'   , 'rebin':3, 'split_prompt_ph':True, 'unblind':True, 'logy':False}
        })
    # for name in ['kin', 'loose', 'medium', 'tight']:
    #     VarInfo_VVGamma.update({
    #         'sigmaiEtaiEta_'+name+'Photons': ['#sigma_{i#etai#eta}', 1, True]
    #     })
    VarInfo_VVGamma.update({
        'kinPh_central_N'       : {'title':'Number of #gamma_{kin}'    , 'split_prompt_ph':False, 'unblind':True , 'logy':True, 'text':False},
        'veryLoosePh_central_N' : {'title':'Number of #gamma_{loose}'  , 'split_prompt_ph':True, 'unblind':False, 'logy':True, 'text':True},
        'loosePh_central_N'     : {'title':'Number of #gamma_{tight}'  , 'split_prompt_ph':True, 'unblind':False, 'logy':True, 'text':True},
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

    VarInfo_VVGamma.update({
        'MET'         : {'title': 'MET [GeV]', 'xmax': 120.},
        'lead_lep_pt' : {},
        'lead_lep_eta': {},
        'GEN_chLeptons' : {'title':'# GEN charged leptons' ,'unblind':False, 'logy':True, 'ymin':1}
    })

    return VarInfo_VVGamma
