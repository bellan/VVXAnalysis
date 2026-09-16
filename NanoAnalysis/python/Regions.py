from enum import IntFlag


# IntFlag is an enum that supports bitwise operations like OR e AND
class Flags(IntFlag):
    ######## Pure flags ########
    
    #### Photons bits 0-7 ####
    P1K       = 1 << 0 # at least one photon kinem. acc. (kin photon) 

    #### Photons cutBased ID wp bits 1-4 ####
    P1cutVL   = 1 << 1 # at least one kin photon passes SMP-24-014 cutBased ID veryLoose wp  
    P1cutL    = 1 << 2 # at least one kin photon passes EGM cutBased ID Loose wp
    P1cutM    = 1 << 3 # at least one kin photon passes EGM cutBased ID Medium wp
    P1cutT    = 1 << 4 # at least one kin photon passes EGM cutBased ID Tight wp 

    #### Photons MVA ID wp bits 5-6 ####
    P1mvaL    = 1 << 5 # at least one kin photon passes EGM MVA ID Loose wp (wp90)
    P1mvaT    = 1 << 6 # at least one kin photon passes EGM MVA ID Tight wp (wp80)

    #### Jets bits 8-15 ####
    J2        = 1 << 8   # 2 or more jets in kinem. acc. pass the JME ID Loose wp
    J3        = 1 << 9   # 3 or more jets in kinem. acc. pass the JME ID Loose wp
    J4        = 1 << 10  # 4 or more jets in kinem. acc. pass the JME ID Loose wp
    J5        = 1 << 11  # 5 or more jets in kinem. acc. pass the JME ID Loose wp

    #### Leptons bits 16-39 ####

    ## 2 leptons bits 16-23 ##
    L2P    = 1 << 19
    L1P1F  = 1 << 18
    L2F    = 1 << 17
    L2L    = 1 << 16

    ## 3 leptons bits 24-31 ##
    L3P    = 1 << 29
    L2P1L  = 1 << 28
    L2P1F  = 1 << 27
    L1P2F  = 1 << 26
    L3F    = 1 << 25
    L3L    = 1 << 24

    ## 4 leptons bits 32-39 ##
    L4P    = 1 << 37
    L3P1F  = 1 << 36
    L2P2F  = 1 << 35
    L1P3F  = 1 << 34
    L4F    = 1 << 33
    L4L    = 1 << 32


    ## Composite flags ##
    # NB: they are *masks*, therefore to compose them we should use "|" instead of "&"

    ## Flags for gamma categorization ##
    P1cutVLfailL       = P1cutVL  | ~(P1cutL )
    P1mvaKfailL        = P1K      | ~(P1mvaL )

    ## Flags for nJet categorization ##
    J2cat              = J2       | ~(J3)
    J3cat              = J3       | ~(J4)
    J4cat              = J4       | ~(J5)
     
    ## Flags for CRL4P_P1F ##
    L4P_P1cutVLfailL   = L4P      | P1cutVLfailL

    ## Flags for SRL4P_P1P ##    
    L4P_P1cutL         = L4P      | P1cutL
    L4P_P1mvaL         = L4P      | P1mvaL

    ## Flags for CRL3P_P1F ##
    L3P_P1cutVLfailL   = L3P      | P1cutVLfailL
    L3P_P1mvaKfailL    = L3P      | P1mvaKfailL

    ## Flags for SRL3P_P1P ##    
    L3P_P1cutL         = L3P      | P1cutL
    L3P_P1mvaL         = L3P      | P1mvaL

    ## Flags for CRL2P_P1F ##
    L2P_P1cutVLfailL   = L2P      | P1cutVLfailL
    L2P_P1mvaKfailL    = L2P      | P1mvaKfailL
    
    ## Flags for SRL2P_P1P ##    
    L2P_P1cutL         = L2P      | P1cutL
    L2P_P1mvaL         = L2P      | P1mvaL

    ## Flags for SRL2P_J2+ ##    
    L2P_J2             = L2P      | J2

    ## Flags for SRL2P_P1P_J2+ ##    
    L2P_P1cutL_J2      = L2P      | P1cutL | J2
    L2P_P1mvaL_J2      = L2P      | P1mvaL | J2

    @staticmethod
    def check(regionWord, region):
        return regionWord & region == region

pt10 = {'name': 'pt10',
        'cuts': {'pt': ('>', 10)},
        'min_particles': 2,
        'max_particles': float('inf')}

pt20 = {'name': 'pt20',
        'cuts': {'pt': ('>', 20)},
        'min_particles': 1,
        'max_particles': float('inf')}

eta4p7 =  {'name': 'eta4p7', 
           'cuts': {'eta' : ('<', 4.7, '||')},
           'min_particles': 1,
           'max_particles': float('inf')}


# TRANSITION_BARREL_ENDCAP = 1.479
# region not covered : [1.4442,1.566]

EBmaxEta =  {'name': 'inEB', 
             'cuts': {'eta' : ('<', 1.4442, '||')},
             'min_particles': 1,
             'max_particles': float('inf')}

EEminEta =  {'name': 'inEE',
             'cuts': {'eta' : ('>', 1.566, '||')},
             'min_particles': 1,
             'max_particles': float('inf')}

EEmaxEta =  {'name': 'EEmaxEta', 
             'cuts': {'eta' : ('<', 2.5, '||')},
             'min_particles': 1,
             'max_particles': float('inf')}


P1K =  {'name': 'KinPhotons',
        'cuts': {'pt' : ('>', 20), 'eta' : ('<', 2.5, '||')},
        'min_particles': 1,
        'max_particles': float('inf')}

P1mvaL = {**P1K, 'name': 'MvaIdLoosePhotons', 'mvaID_WP90' : ('==', True)}

"""
P1cutVL =  {'name': 'CutIdVeryLoosePhotons',
           'cuts': {'pt' : ('>', 20), 'eta' : ('<', 2.5, '||'), 'hoe' : ('>', 20)},
           'min_particles': 1,
           'max_particles': float('inf')}
"""
P1cutVL = {**P1K, 'name': 'CutIdVeryLoosePhotons'} # FIXME! Currently left as P1K
P1cutL  = {**P1K, 'name': 'CutIdLoosePhotons',  'cutBased' : ('>=', 1)}
P1cutM  = {**P1K, 'name': 'CutIdMediumPhotons', 'cutBased' : ('>=', 2)}
P1cutT  = {**P1K, 'name': 'CutIdTightPhotons',  'cutBased' : ('>=', 3)}

J2 = {'name': 'J2',
      'cuts': {'pt' : ('>', 20), 'eta' : ('<', 4.7, '||'), 'jetId' : ('>=', 2)},
      'min_particles': 2,
      'max_particles': float('inf')}

J3    = {**J2, 'name': 'J3', 'min_particles': 3}
J4    = {**J2, 'name': 'J4', 'min_particles': 4}
J5    = {**J2, 'name': 'J5', 'min_particles': 5}
J2cat = {**J2, 'name': 'J2cat', 'max_particles': 2}
J3cat = {**J3, 'name': 'J3cat', 'max_particles': 3}
J4cat = {**J4, 'name': 'J4cat', 'max_particles': 4}

L1P = {'name': 'TightLeptons',
       'cuts': {'ZZFullSel' : ('==', True)},
       'min_particles': 1,
       'max_particles': 1}

L1L = {'name': 'LooseLeptons',
       'cuts': {'ZZRelaxedId' : ('==', True)},
       'min_particles': 1,
       'max_particles': 1}

L1F = {'name': 'FailLeptons',
       'cuts': {'ZZFullSel' : ('==', False)},
       'min_particles': 1,
       'max_particles': 1}

L4P = {**L1P, 'name': 'L4P', 'min_particles': 4, 'max_particles': 4}
L3P = {**L1P, 'name': 'L3P', 'min_particles': 3, 'max_particles': 3}
L2P = {**L1P, 'name': 'L2P', 'min_particles': 2, 'max_particles': 2}


L4F = {**L1F, 'name': 'L4F', 'min_particles': 4, 'max_particles': 4}
L3F = {**L1F, 'name': 'L3F', 'min_particles': 3, 'max_particles': 3}
L2F = {**L1F, 'name': 'L2F', 'min_particles': 2, 'max_particles': 2}

L4L = {**L1L, 'name': 'L4L', 'min_particles': 4, 'max_particles': 4}
L3L = {**L1L, 'name': 'L3L', 'min_particles': 3, 'max_particles': 3}
L2L = {**L1L, 'name': 'L2L', 'min_particles': 2, 'max_particles': 2}

# Special, for ZL region --> L2P1L
L23P = {**L1P, 'name': 'L23P', 'min_particles': 2, 'max_particles': 3} 


flagDefinitions = [
    ########################################################
    {'name'   : 'L4P',        
     'leptons'  : {'selection': [pt10, pt20, L4P]}
     },
    ########################################################

    ########################################################
    {'name'   : 'L3P1F',        
     'leptons'  : {'selection': [pt10, pt20, L4L, L3P, L1F]}
     },
    
    ########################################################
    {'name'   : 'L2P2F',        
     'leptons'  : {'selection': [pt10, pt20, L4L, L2P, L2F]}
     },

    ########################################################
    {'name'   : 'L4F',        
     'leptons'  : {'selection': [pt10, pt20, L4L, L4F]}
     },
    ########################################################

    ########################################################
    {'name'   : 'L4L',        
     'leptons'  : {'selection': [pt10, pt20, L4L]}
     },

    ########################################################
    {'name'   : 'L3P',        
     'leptons'  : {'selection': [pt10, pt20, L3P]}
     },

    ########################################################
    {'name'   : 'L2P1L',        
     'leptons'  : {'selection': [pt10, pt20, L3L, L23P]}
     },
    
    ########################################################
    {'name'   : 'L2P1F',        
     'leptons'  : {'selection': [pt10, pt20, L3L, L2P, L1F]}
     },

    ########################################################
    {'name'   : 'L1P2F',        
     'leptons'  : {'selection': [pt10, pt20, L3L, L1P, L2F]}
     },

    ########################################################
    {'name'   : 'L3F',        
     'leptons'  : {'selection': [pt10, pt20, L3L, L3F]}
     },
    
    ########################################################
    {'name'   : 'L3L',        
     'leptons'  : {'selection': [pt10, pt20, L3L]}
     },

    ########################################################
    {'name'   : 'L2P',        
     'leptons'  : {'selection': [pt10, pt20, L2P]}
     },

    ########################################################
    {'name'   : 'L1P1F',        
     'leptons'  : {'selection': [pt10, pt20, L2L, L1P, L1F]}
     },

    ########################################################
    {'name'   : 'L2F',        
     'leptons'  : {'selection': [pt10, pt20, L2L, L2F]}
     },

    ########################################################
    {'name'   : 'L2L',        
     'leptons'  : {'selection': [pt10, pt20, L2L]}
     },
    
    ########################################################
    {'name'   : 'P1cutL',        
     'photons'  : {'selection': [pt20, EEmaxEta, P1cutL]}
     },
    
    ########################################################
    {'name'   : 'P1mvaL',        
     'photons'  : {'selection': [pt20, EEmaxEta, P1mvaL]}
     },
   
    ########################################################    
    {'name' : 'J2',
     'jets'    : {'selection': [pt20, eta4p7, J2]}
     },

]

