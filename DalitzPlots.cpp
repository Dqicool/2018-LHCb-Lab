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

    #define DRAW_BAC
    #define DRAW_COM
    #define DRAW_SIG
    #define DRAW_ASY

// Define Data Choice
    #define DATA_ALL
    //#define DATA_UP
    //#define DATA_DOWN
    //#define DATA_PHASE_SPACE
//Main function
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
        //7.3 getting signal events positive
            TH2F *DalitzPosCom = (TH2F*) f->Get("h_Dalitz_Pos_Com");
            TH2F *DalitzPosBac = (TH2F*) f->Get("h_Dalitz_Pos_Bac");
        //7.4 getting signal events negtive
            TH2F *DalitzNegCom = (TH2F*) f->Get("h_Dalitz_Neg_Com");
            TH2F *DalitzNegBac = (TH2F*) f->Get("h_Dalitz_Neg_Bac");

        //rescale 
            for(int i = 0; i<100; i++)
            {
                for( int j=0;j<100;j++)
                {
                    DalitzNegBac->SetBinContent(i,j,DalitzNegBac->GetBinContent(i,j)/3);
                    DalitzPosBac->SetBinContent(i,j,DalitzPosBac->GetBinContent(i,j)/3);
                }
            } 
        //extract signal
            TH2F *DalitzPosSig = (TH2F*) DalitzPosCom->Clone("DalitzPosSig");
            DalitzPosSig->Add(DalitzPosBac,-1);
            TH2F *DalitzNegSig = (TH2F*) DalitzNegCom->Clone("DalitzNegSig");
            DalitzNegSig->Add(DalitzNegBac,-1);

        //7.5 Calculate Global asymmytries
            TH2F *NumAll = (TH2F*) DalitzPosSig->Clone("NumAll");
            NumAll->Add(DalitzNegSig);
            Double_t NumB = 0;
            for(int i = 0; i<100; i++)
            {
                for( int j=0;j<100;j++)
                {
                    NumB += NumAll->GetBinContent(i,j);
                }
            } 
            cout<<NumB<<endl;

            TH2F *NumMinus = (TH2F*)DalitzNegSig->Clone("NumMinus");
            NumMinus->Add(DalitzPosSig,-1);
            Double_t NumA=0;
            for(int i = 0; i<100; i++)
            {
                for( int j=0;j<100;j++)
                {
                    NumA += NumMinus->GetBinContent(i,j);
                }
            } 
            cout<<NumA<<endl;
            cout<<NumA/NumB<<endl;
            //Calculating Local Asymmetry
            TH2F *ASSY = (TH2F*)DalitzNegSig->GetAsymmetry(DalitzPosSig);
        //Drawing
            #ifdef DRAW_COM
                TCanvas *c10 = new TCanvas("c10","",CANVASIZE,CANVASIZE);
                DalitzPosCom->SetAxisRange(0,500,"Z");
                DalitzPosCom->Draw("colz");
                DalitzPosCom->Draw("TEXT SAME");
                c10->SaveAs("Plots/c10_DalitzPlot_Pos_Com.pdf");

                TCanvas *c13 = new TCanvas("c13","",CANVASIZE,CANVASIZE);
                DalitzNegCom->SetAxisRange(0,500,"Z");
                DalitzNegCom->Draw("colz");
                DalitzNegCom->Draw("TEXT SAME");
                c13->SaveAs("Plots/c13_DalitzPlot_Neg_Com.pdf");
            #endif

            #ifdef DRAW_BAC
                TCanvas *c14 = new TCanvas("c14","",CANVASIZE,CANVASIZE);
                DalitzNegBac->SetAxisRange(0,500,"Z");
                DalitzNegBac->Draw("colz");
                DalitzNegBac->Draw("TEXT SAME");
                c14->SaveAs("Plots/c14_DalitzPlot_Neg_Bac.pdf");

                TCanvas *c11 = new TCanvas("c11","",CANVASIZE,CANVASIZE);
                DalitzPosBac->SetAxisRange(0,500,"Z");
                DalitzPosBac->Draw("colz");
                DalitzPosBac->Draw("TEXT SAME");
                c11->SaveAs("Plots/c11_DalitzPlot_Pos_Bac.pdf");
            #endif


            #ifdef DRAW_SIG
                TCanvas *c12 = new TCanvas("c12","",CANVASIZE,CANVASIZE);
                DalitzPosSig->SetAxisRange(0,500,"Z");
                DalitzPosSig->Draw("colz");
                DalitzPosSig->Draw("TEXT SAME");
                c12->SaveAs("Plots/c12_DalitzPlot_Pos_Sig.pdf");

                TCanvas *c15 = new TCanvas("c15","",CANVASIZE,CANVASIZE);
                DalitzNegSig->SetAxisRange(0,500,"Z");
                DalitzNegSig->Draw("colz");
                DalitzNegSig->Draw("TEXT SAME");
                c15->SaveAs("Plots/c15_DalitzPlot_Pos_Sig.pdf");
            #endif

            #ifdef DRAW_ASY

                TCanvas *c16 = new TCanvas("c16","",CANVASIZE,CANVASIZE);
                ASSY->SetAxisRange(-1,1,"Z");
                ASSY->Draw("colz1");
                ASSY->Draw("TEXT SAME");
                ASSY->SetBarOffset(0.2);
                NumMinus->SetBarOffset(-0.2);
                NumMinus->Draw("SAME TEXT ");
                c16->SaveAs("Plots/c16_Local_Assymmytry.pdf");

                TCanvas *c17 = new TCanvas("c17","",CANVASIZE,CANVASIZE);
                ASSY->SetAxisRange(-1,1,"Z");
                NumMinus->Draw("colz1");
                ASSY->Draw("TEXT SAME");
                ASSY->SetBarOffset(0.2);
                NumMinus->SetBarOffset(-0.2);
                NumMinus->Draw("SAME TEXT ");
                c17->SaveAs("Plots/c17_Local_ABS_Assymmytry.pdf");
            #endif
    }