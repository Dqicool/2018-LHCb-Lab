//1 include 
    #include <TFile.h>
    #include <TH1.h>
    #include <TH2.h>
    #include <TCanvas.h>
    #include <TMath.h>
    #include <TF1.h>
    #include <TF2.h>
    #include <iostream>
    using namespace std;
// Define Data Choice
    #define DATA_ALL
    //#define DATA_UP
    //#define DATA_DOWN
    //#define DATA_PHASE_SPACE

void DalitzPlots(){
        //7.2 declare Data choice
            #ifdef DATA_ALL
                cout<<"Using all Data"<<endl;
            #endif
            #ifdef DATA_UP
                cout<<"Using Magnet Up Data"<<endl;
            #endif
            #ifdef DATA_DOWN
                cout<<"Using Magnet Down Data"<<endl;
            #endif
            #ifdef DATA_PHASE_SPACE
                cout<<"Using Phase Space Data"<<endl;
            #endif
        //7.3 load chosen Data from file
            #ifdef DATA_ALL
                TFile *f = new TFile("Output/DataAll.root");
		        TH2F *DalitzCom = (TH2F*) f->Get("h_Dalitz_Com");
                TH2F *DalitzBac = (TH2F*) f->Get("h_Dalitz_Bac");
                TH2F *DalitzSig = (TH2F*) DalitzCom->Clone("DalitSig");
                DalitzSig->Add(DalitzBac,-1);
            #endif
            #ifdef DATA_DOWN
                TFile *f = new TFile("Output/DataMagnetDown.root");
		        TH2F *Dalitz = (TH2F*) f->Get("h_Dalitz");
            #endif
            #ifdef DATA_UP
                TFile *f = new TFile("Output/DataMagnetUp.root");
		        TH2F *Dalitz = (TH2F*) f->Get("h_Dalitz");
            #endif
            #ifdef DATA_PHASE_SPACE
                TFile *f = new TFile("Output/PhaseSpace.root");
		        TH2F *Dalitz = (TH2F*) f->Get("h_Dalitz");
            #endif
//7.3 

TCanvas *c10 = new TCanvas("c10","",600,600);
DalitzCom->Draw("colz");
c10->SaveAs("Plots/c10_DalitzPlot_Com.pdf");

TCanvas *c11 = new TCanvas("c11","",600,600);
DalitzBac->Draw("colz");
c11->SaveAs("Plots/c11_DalitzPlot_Bac.pdf");

TCanvas *c12 = new TCanvas("c12","",600,600);
DalitzSig->Draw("colz");
c12->SaveAs("Plots/c12_DalitzPlot_Sig.pdf");

}