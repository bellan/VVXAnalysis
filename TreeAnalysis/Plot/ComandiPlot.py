#Comandi per lanciare Stack_func.py per i vari plot richesti


#Plot 1.
# SR: bkg irriducibile tutto insieme
#     bkg riducibile 
#     segnale diviso nelle 3 componenti (2 stato iniziale + other)

./scriptStack.py SR sig+irred_bkg_tot

#Plot 2.
# SR: bkg irriducibile diviso nelle 3 componenti

./scriptStack.py SR irred_bkg_div (commentare fake rate nello stack plot, setmaximum, data)

#Plot 3.
# SR: bkg riducibile da MC diviso per i vari campioni
#     fake rate sovrapposto

./scriptStack.py SR_compare red_bkg (setmaximum, data)

#Plot 4.
# SR: bkg riducibile 
#     bkg irriducibile tutto insieme
#     segnale diviso per stati finali

./scriptStack.py SR sig+irred_bkg_tot (il irred_bkg è già incluso, cambiare funzione, non setmaximum)
NB: CAMBIARE IL NOME AL FILE SALVATO SE NO SI SOVRAPPONE!!!

#Plot 5.
# SR: campioni di segnale (tutto oppure solo i due maggiori contributi qq,gg) diviso in due parti, se passa la definizione di segnale o no

./scriptStack.py SR sig (sig_qq, sig_gg) (cambiare funzione)
NB: MEGLIO CAMBIARE IL NOME AL FILE SALVATO 
