from enum import IntFlag

# IntFlag is an enum that supports bitwise operations like OR e AND
class Regions(IntFlag):
    R4P    = 1 << 0  
#    R3P1F  = 1 << 1  
#    R2P1F  = 1 << 2
#    R4P_1P
#    R4P_1L
#    R4P_1F
    
    
    R3P   = 1 << 1  # Valore: 2  (bit 1)
#    R2P1F
#    R1P2F
#    R3F

    R2P   = 1 << 2  # Valore: 4  (bit 2)

    

 #   R2P1L
