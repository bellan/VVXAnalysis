/**
 *  Macro to get efficiency of Bosons reconstruction from results of VZGAnalyzer 
 *
 *
 *	Usage: 	root [-l] [-b] [-q] 'VZGammaRecoEfficiencyAnalysis.C("<sample name>")'		-l do not display banner 	-b run in background		-q close after finishing
 *	e.g.: 	root -l 'VZGammaRecoEfficiencyAnalysis.C("WZGTo2L2jG")'
 *	It may become necessary to change the path and/or to expand the samples list	
 *			
 *  $Date: 2023/06/26 10:17:08 
 *  $Revision: 2.0 $
 *
 *  \author C. Tarricone cristiano.tarricone@cern.ch
 */

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <iostream>
#include <vector>

TString resDir = "resultsForReweighting";
const std::vector<TString> years = {"2016preVFP", "2016postVFP", "2017", "2018"};

void DEBUGextractWeights() {
  gStyle->SetOptStat(0);

  // Vettore di file fuori dal loop
    std::vector<TFile*> DYFiles;
    std::vector<TFile*> DataFiles;
    std::vector<TFile*> ZGFiles;
    std::vector<TFile*> TTFiles;

    Color_t colors[] = {kRed+1, kSpring+3, kBlue, kViolet  };
    Int_t markerStyles[] = {20,22,23,24};
    TLegend* legend = new TLegend(0.12, 0.68, 0.32, 0.88);

    // Ciclo per aprire i file prima del loop
    for (int i = 0; i < 4; i++) {
        TString year = years[i];
        TString path = "../" + resDir + "/" + year + "/VZGAnalyzer_SR2P/";

        TFile* DYFile = new TFile(path + "DYJetsToLL_M50.root", "READ");
        TFile* DataFile = new TFile(path + "data_obs.root", "READ");
        TFile* ZGFile = new TFile(path + "ZGToLLG.root", "READ");
        TFile* TTFile = new TFile(path + "TTTo2L2Nu.root", "READ");

        if (!DYFile->IsOpen() || !DataFile->IsOpen()) {
            std::cerr << "Errore nell'aprire i file: " << path << std::endl;
            continue;
        }

        // Aggiungi i file ai rispettivi vettori
        DYFiles.push_back(DYFile);
        DataFiles.push_back(DataFile);
        ZGFiles.push_back(ZGFile);
        TTFiles.push_back(TTFile);
    }

    // Ora che i file sono aperti, iniziamo a processarli
    std::vector<TH1F*> hEffs;  // Vettore per memorizzare gli istogrammi

    for (int i = 0; i < 4; i++) {
        TString year = years[i];
        std::cout << "--------------- " << year << " ---------------" << std::endl;

        // Recupera i file dal vettore
        TFile* DYFile = DYFiles[i];
        TFile* DataFile = DataFiles[i];
        TFile* ZGFile = ZGFiles[i];
        TFile* TTFile = TTFiles[i];

        TH1F* hNum = (TH1F*)DataFile->Get("REWGT_numDATA");
        TH1F* hDen = (TH1F*)DYFile->Get("REWGT_DY-GEN");
        TH1F* hZG = (TH1F*)ZGFile->Get("REWGT_else-GEN");
        TH1F* hTT = (TH1F*)TTFile->Get("REWGT_else-GEN");
	hNum->Add(hZG, -1);
	hNum->Add(hTT, -1);
	cout <<"{ ";
	for (int j = 1; j <= hDen->GetNbinsX(); j++) {
	  cout << hNum->GetBinContent(j)/hDen->GetBinContent(j)<<", ";

	  if (hDen->GetBinContent(j) <= 0) {
	    cout << "Warning: Bin " << j << " in hDen is zero!" << endl;
	    hDen->SetBinContent(j, 1e-6); // Evita divisioni per zero
	  }
	}
	cout <<"} "<<endl;

	/*
        if (!hNum || !hDen) {
            std::cerr << "Errore nella lettura degli istogrammi." << std::endl;
            continue;
        }

        std::cout << "Numero di entry in hNum: " << hNum->GetEntries() << std::endl;
        std::cout << "Numero di entry in hDen: " << hDen->GetEntries() << std::endl;
        hNum->Sumw2();
        hDen->Sumw2();

        // Verifica che i dati siano validi prima della divisione
        if (hDen->GetEntries() == 0) {
            std::cerr << "Errore: Denominatore hDen ha 0 entry!" << std::endl;
	continue;
        }
	*/
        // Clona hNum per creare hEff, separato dalla memoria legata al file
        TH1F* hEff = (TH1F*)hNum->Clone("hEff"); // Cloniamo hNum in un nuovo oggetto hEff

        hEff->Divide(hDen);
	/*
        // Verifica se ci sono valori invalidi nell'istogramma risultante
        for (int bin = 1; bin <= hEff->GetNbinsX(); bin++) {
            if (std::isnan(hEff->GetBinContent(bin))) {
                std::cerr << "Attenzione: Valore NaN in hEff, bin " << bin << std::endl;
            }
        }
	
        std::cout << "Numero di entry in hEff dopo la divisione: " << hEff->GetEntries() << std::endl;
	*/
        hEffs.push_back(hEff);  // Aggiungi l'oggetto hEff al vettore
    }

    // Verifica che il vettore hEffs contenga oggetti validi
    if (hEffs.empty()) {
        std::cerr << "hEffs è vuoto! Controlla che gli istogrammi siano stati caricati correttamente." << std::endl;
        return;
    }

    // Canvas per visualizzare solo un istogramma
    TCanvas* cDrawing = new TCanvas("weights vs Z pt ", "weights vs Z pt ", 0, 0, 800, 800);

    // Disegna solo il primo istogramma per isolare il problema
    /*
      if (hEffs[0]) {
        std::cout << "Disegnando il primo istogramma hEffs[0]..." << std::endl;
        hEffs[0]->SetLineColor(kRed);
        hEffs[0]->SetMarkerStyle(20);
        hEffs[0]->SetMarkerColor(kRed);

        hEffs[0]->Draw("AP");

        hEffs[0]->GetXaxis()->SetTitle("p_{T} Z Cand. [GeV]");
        hEffs[0]->GetXaxis()->SetTitleSize(0.035);
        hEffs[0]->GetYaxis()->SetTitle("weights");
        hEffs[0]->GetYaxis()->SetTitleSize(0.035);
    } else {
        std::cerr << "Errore: hEffs[0] è nullo!" << std::endl;
    }
    */
    for(int i = 0; i<4; i++){
      hEffs[i]->SetLineColor(colors[i]);
      hEffs[i]->SetMarkerStyle(markerStyles[i]);
      hEffs[i]->SetMarkerColor(colors[i]);

      hEffs[i]->Draw(i == 0 ? "E" : "E SAME");
      if(i==0){
	hEffs[i]->GetXaxis()->SetLabelSize(0.035); // Imposta la dimensione delle etichette per l'asse X
	hEffs[i]->GetYaxis()->SetLabelSize(0.035); // Imposta la dimensione delle etichette per l'asse Y
	hEffs[i]->GetXaxis()->SetTitleSize(0.035); // Imposta la dimensione del titolo per l'asse X
	hEffs[i]->GetYaxis()->SetTitleSize(0.035); // Imposta la dimensione del titolo per l'asse Y
      }

      hEffs[i]->GetXaxis()->SetTitle("p_{T} Z Cand. [GeV]");
      hEffs[i]->GetXaxis()->SetTitleSize(0.035);
      hEffs[i]->GetYaxis()->SetTitle("weights");
      hEffs[i]->GetYaxis()->SetTitleSize(0.035);
      if (legend) {
        legend->SetBorderSize(1);
        legend->AddEntry(hEffs[i], years[i] + " weights", "lp");
        legend->Draw();
      } else {
        std::cerr << "Errore nell'inizializzare la legenda." << std::endl;
      }

    }
      
    // Verifica che la legenda è correttamente inizializzata
      
    // Verifica che il canvas sia correttamente modificato e aggiornato
    cDrawing->Modified();
    cDrawing->Update();

    TPaveText *paveTLeft = new TPaveText(0.065, 0.87, 0.9, 0.95, "NDCNDC");
    paveTLeft->SetFillColor(0);
    paveTLeft->SetFillStyle(0);
    paveTLeft->SetBorderSize(0);
    paveTLeft->SetTextAlign(11);
    paveTLeft->SetTextFont(62);
    paveTLeft->SetTextSize(0.04);
    paveTLeft->AddText("CMS #bf{#it{Internal}}");

    TPaveText *paveTRight = new TPaveText(0.77, 0.87, 0.87, 0.95, "NDCNDC");
    paveTRight->SetFillColor(0);
    paveTRight->SetFillStyle(0);
    paveTRight->SetBorderSize(0);
    paveTRight->SetTextAlign(11);
    paveTRight->SetTextFont(62);
    paveTRight->SetTextSize(0.035);
    paveTRight->AddText("#bf{13 TeV}");

    cDrawing->cd();
    paveTLeft->Draw();
    paveTLeft->SetBit(kCanDelete);
    paveTRight->Draw();
    paveTRight->SetBit(kCanDelete);

    // Salva l'immagine
    TString fname = "Zpt_reweight";
    cDrawing->SaveAs("efficiencyPlots/" + fname + ".png");
    cDrawing->SaveAs("efficiencyPlots/" + fname + ".pdf");

    std::cout << "Operazione completata senza crash." << std::endl;

    // Liberazione della memoria
    for (TH1F* h : hEffs) {
        delete h;  // Assicurati di liberare la memoria degli oggetti TH1F
    }

        // Chiudi i file dopo aver finito di elaborarli
    for (TFile* DYFile : DYFiles) {
        DYFile->Close();
        delete DYFile;
    }
    for (TFile* DataFile : DataFiles) {
        DataFile->Close();
        delete DataFile;
    }
    for (TFile* ZGFile : ZGFiles) {
        ZGFile->Close();
        delete ZGFile;
    }
    for (TFile* TTFile : TTFiles) {
        TTFile->Close();
        delete TTFile;
    }

}
