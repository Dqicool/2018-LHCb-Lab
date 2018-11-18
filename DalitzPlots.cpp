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
    #define CANVASIZE 400
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
            #endif
            #ifdef DATA_DOWN
                TFile *f = new TFile("Output/DataMagnetDown.root");
            #endif
            #ifdef DATA_UP
                TFile *f = new TFile("Output/DataMagnetUp.root");
            #endif
            #ifdef DATA_PHASE_SPACE
                TFile *f = new TFile("Output/PhaseSpace.root");
            #endif
//7.3 
    TH2F *DalitzPosCom = (TH2F*) f->Get("h_Dalitz_Pos_Com");
    TH2F *DalitzPosBac = (TH2F*) f->Get("h_Dalitz_Pos_Bac");
    TH2F *DalitzPosSig = (TH2F*) DalitzPosCom->Clone("DalitPosSig");
    DalitzPosSig->Add(DalitzPosBac,-1);
    /*//!Set All Minus Bins To 0 
    for(int i = 0; i<100; i++){
        for( int j=0;j<100;j++){
            DalitzPosSig->SetBinContent(i,j,DalitzPosSig->GetBinContent(i,j)<0?0:DalitzPosSig->GetBinContent(i,j));
        }
    } */

    TH2F *DalitzNegCom = (TH2F*) f->Get("h_Dalitz_Neg_Com");
    TH2F *DalitzNegBac = (TH2F*) f->Get("h_Dalitz_Neg_Bac");
    TH2F *DalitzNegSig = (TH2F*) DalitzPosCom->Clone("DalitNegSig");
    DalitzNegSig->Add(DalitzNegBac,-1);
    /*//!Set All Minus Bins To 0 
    for(int i = 0; i<100; i++){
        for( int j=0;j<100;j++){
            DalitzNegSig->SetBinContent(i,j,DalitzNegSig->GetBinContent(i,j)<0?0:DalitzNegSig->GetBinContent(i,j));
        }
    } */

TCanvas *c10 = new TCanvas("c10","",CANVASIZE,CANVASIZE);
DalitzPosCom->SetAxisRange(0,500,"Z");
DalitzPosCom->Draw("colz");
c10->SaveAs("Plots/c10_DalitzPlot_Pos_Com.pdf");

TCanvas *c11 = new TCanvas("c11","",CANVASIZE,CANVASIZE);
DalitzPosBac->SetAxisRange(0,500,"Z");
DalitzPosBac->Draw("colz");
c11->SaveAs("Plots/c11_DalitzPlot_Pos_Bac.pdf");

TCanvas *c12 = new TCanvas("c12","",CANVASIZE,CANVASIZE);
DalitzPosSig->SetAxisRange(0,500,"Z");
DalitzPosSig->Draw("colz");
c12->SaveAs("Plots/c12_DalitzPlot_Pos_Sig.pdf");

TCanvas *c13 = new TCanvas("c13","",CANVASIZE,CANVASIZE);
DalitzNegCom->SetAxisRange(0,500,"Z");
DalitzNegCom->Draw("colz");
c13->SaveAs("Plots/c13_DalitzPlot_Pos_Com.pdf");

TCanvas *c14 = new TCanvas("c14","",CANVASIZE,CANVASIZE);
DalitzNegBac->SetAxisRange(0,500,"Z");
DalitzNegBac->Draw("colz");
c14->SaveAs("Plots/c14_DalitzPlot_Pos_Bac.pdf");

TCanvas *c15 = new TCanvas("c15","",CANVASIZE,CANVASIZE);
DalitzNegSig->SetAxisRange(0,500,"Z");
DalitzNegSig->Draw("colz");
c15->SaveAs("Plots/c15_DalitzPlot_Pos_Sig.pdf");

}