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
		TH2F *Dalitz = (TH2F*) f->Get("h_Dalitz");
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
Dalitz->Draw("colz");
c10->SaveAs("Plots/c10_DalitzPlot.pdf");

}