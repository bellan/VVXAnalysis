from enum import IntFlag


# IntFlag is an enum that supports bitwise operations like OR e AND
class Flags(IntFlag):
    ######## Pure flags ########

    #### Photons bits 0-7 ####
    P1P    = 1 << 2
    P1F    = 1 << 1
    P1L    = 1 << 0

    #### Jets bits 8-15 ####
    J2    = 1 << 8

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
    L4P_L1P  = L4P | P1P
    L4P_J2   = L4P | J2
    L3P_P1P  = L3P | P1P
    L2P_J2   = L2P | J2
    L2P_P1P  = R2P | P1P

    
pt10 = {'name': 'pt10',
        'cuts': {'pt': ('>', 10)},
        'min_particles': 2,
        'max_particles': float('inf')}

pt20 = {'name': 'pt20',
        'cuts': {'pt': ('>', 20)},
        'min_particles': 1,
        'max_particles': float('inf')}

P1P =  {'name': 'TightPhotons',
        'cuts': {'mvaID' : ('>', 0)}, #FIXME!!!
        'min_particles': 1,
        'max_particles': float('inf')}

J2 = {'name': 'J2',
         'cuts': {'pt' : ('>', 20), 'eta' : ('<', 4.7, '||')},
         'min_particles': 2,
         'max_particles': float('inf')}


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
    {'name'   : 'P1P',        
     'photons'  : {'selection': [pt20, P1P]}
     },

   
    
    ########################################################
    # {'name'   : 'P1F',        
    #  'photons'  : {'selection': [pt20, P1L, P1F]}
    #  },

    # ########################################################
    # {'name'   : 'P1L',        
    #  'photons'  : {'selection': [pt20, P1P]}
    #  },
    
    ########################################################
    {'name' : 'J2',
     'jets'    : {'selection': [J2]}
     },

    
]
