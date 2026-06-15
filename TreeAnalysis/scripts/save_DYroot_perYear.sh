#!/bin/bash

# Controlla se l'utente ha fornito l'argomento (nome della directory)
if [ $# -ne 1 ]; then
    echo "Uso: $0 <nome_cartella_output>"
    exit 1
fi

# Prende il primo argomento come nome della cartella
output_dir="$1"

# Percorso del file di output
2016preVFP_output_file = "$output_dir/2016preVFP/VZGAnalyzer_SR2P/DYJetsToLL_M50.root"
2016postVFP_output_file= "$output_dir/2016postVFP/VZGAnalyzer_SR2P/DYJetsToLL_M50.root"
2017_output_file = "$output_dir/2017/VZGAnalyzer_SR2P/DYJetsToLL_M50.root"
2018_output_file = "$output_dir/2018/VZGAnalyzer_SR2P/DYJetsToLL_M50.root"

# Elenco dei file ROOT da unire
2016preVFP_root_files=(
    "results/2016preVFP/2016preVFP/VZGAnalyzer_SR2P/DYJetsToLL_M50_part1of6.root"
    "results/2016preVFP/2016preVFP/VZGAnalyzer_SR2P/DYJetsToLL_M50_part2of6.root"
    "results/2016preVFP/2016preVFP/VZGAnalyzer_SR2P/DYJetsToLL_M50_part3of6.root"
    "results/2016preVFP/2016preVFP/VZGAnalyzer_SR2P/DYJetsToLL_M50_part4of6.root"
    "results/2016preVFP/2016preVFP/VZGAnalyzer_SR2P/DYJetsToLL_M50_part5of6.root"
    "results/2016preVFP/2016preVFP/VZGAnalyzer_SR2P/DYJetsToLL_M50_part6of6.root"
)
2016postVFP_root_files=(
    "results/2016postVFP/2016postVFP/VZGAnalyzer_SR2P/DYJetsToLL_M50_part1of6.root"
    "results/2016postVFP/2016postVFP/VZGAnalyzer_SR2P/DYJetsToLL_M50_part2of6.root"
    "results/2016postVFP/2016postVFP/VZGAnalyzer_SR2P/DYJetsToLL_M50_part3of6.root"
    "results/2016postVFP/2016postVFP/VZGAnalyzer_SR2P/DYJetsToLL_M50_part4of6.root"
    "results/2016postVFP/2016postVFP/VZGAnalyzer_SR2P/DYJetsToLL_M50_part5of6.root"
    "results/2016postVFP/2016postVFP/VZGAnalyzer_SR2P/DYJetsToLL_M50_part6of6.root"
)

# Crea la directory di output se non esiste
if [ ! -d "$output_dir" ]; then
    echo "Creazione della cartella: $output_dir"
    mkdir -p "$output_dir"
fi

# Controlla che tutti i file di input esistano
for file in "${root_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "Errore: il file $file non esiste!"
        exit 1
    fi
done

# Esegue hadd per combinare i file .root
echo "Unendo i file ROOT in $output_file..."
hadd -f "$output_file" "${root_files[@]}"

# Controlla se hadd è andato a buon fine
if [ $? -eq 0 ]; then
    echo "Unione completata con successo! File salvato in: $output_file"
else
    echo "Errore durante l'esecuzione di hadd!"
    exit 1
fi
