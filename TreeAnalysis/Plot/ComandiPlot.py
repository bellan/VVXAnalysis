#Comandi per lanciare scriptStack.py per i vari plot richesti


# ========= Plot 1. ==========

# SR: bkg irriducibile tutto insieme
#     bkg riducibile 
#     segnale diviso nelle 3 componenti (2 stato iniziale + other)

./scriptStack.py SR sig+irred_bkg_tot

# ========== Plot 2. ==========

# SR: bkg riducibile 
#     bkg irriducibile tutto insieme
#     segnale diviso per stati finali

# NB: Cosa cambiare:
#- Cambiare funzione, usare GetStackPlot_fstate nello script
#- Cambiare il nome al file salvato se no si sovrappone!!

./scriptStack.py SR sig+irred_bkg_tot 


# ========= Plot 3. ==========

# SR: bkg irriducibile diviso nelle 3 componenti

# NB: Cosa cambiare:
# - Usare di nuovo la funzione GetStackPlot
# - Commentare "fake rate" nella funzione GetStackPlot
# - Commetare GetDataPlot nllo script
# - Commenatare SetMaximum
# - RICAMBIARE il nome al file salvato!!!!

./scriptStack.py SR irred_bkg_div 

# ========= Plot 4. ==========

# SR: bkg riducibile da MC diviso per i vari campioni
#     fake rate sovrapposto

# NB: Cosa cambiare:
# - Nulla rispetto al Plot 3. 

./scriptStack.py SR_compare red_bkg (setmaximum, data)



# ========= Plot 5.  ==========

# SR: campioni di segnale (tutto oppure solo i due maggiori contributi qq,gg) diviso in due parti, se passa la definizione di segnale o no

# NB: Cosa cambiare:
# - Usare la funzione GetSignalDefPlot

./scriptStack.py SR sig (sig_qq, sig_gg)


# ========= Plot 6,7.  =========

# Control region

# NB: Cosa cambiare:
# - Usare di nuovo la funzione GetStackPlot
# - Rimettere i dati
# - Usare SetMaximum(0.27)

./scriptStack.py CR2P2F sigTot+irred_bkg_tot+red_bkg
./scriptStack.py CR3P1F sigTot+irred_bkg_tot+red_bkg #(per questa regione SetMaximum(0.9)

# ========= PER I PLOT CON FAKE RATE MC DRIVEN  ==========
 
# - Mettere "MC" nel secondo campo ogni volta che GetFakeRate viene chiamato sia in Stack.py che in scriptStack.py
# - Cambiare il nome nella legenda a seconda del metodo in scriptStack.py "Red background from FR method, MC(data) driven"


# ========= PER IL PLOT CONFRONTO FAKE RATE MC/DATA DRIVEN  ==========

# Esiste un altro script (per comodità) chiamato "scriptStack_comparison.py" che deve essere fatto girare solo definendo la "region" come SR_compare
