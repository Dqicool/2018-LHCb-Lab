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
    #define CANVASIZE1 600
    #define CANVASIZE2 400

    #define DRAW_BAC
    #define DRAW_COM
    #define DRAW_SIG
    #define DRAW_ASY

    #define MERGE

// Define Data Choice
    #define DATA_ALL
    //#define DATA_UP
    //#define DATA_DOWN
    //#define DATA_PHASE_SPACE
// Merge Dalitz 
    void MergeBins(TH2F* Histo1, TH2F* Histo2, TH2F* Histo3, TH2F* Histo4,TH2F* Histo5,TH2F* Histo6, double x, double w, double y, double h){
        int xx = 4*x+1; int yy = 4*y+1; int ww = 4*(w-x); int hh = 4*(h-y);
        Double_t sum1 = 0, avg1;
        Double_t sum2 = 0, avg2;
        Double_t sum3 = 0, avg3;
        Double_t sum4 = 0, avg4;
        Double_t sum5 = 0, avg5;
        Double_t sum6 = 0, avg6;
        for (int i = xx; i < xx+ww; i++)
        {
            for (int j = yy; j<yy+hh; j++)
            {
                sum1 += Histo1->GetBinContent(i,j);
                sum2 += Histo2->GetBinContent(i,j);
                sum3 += Histo3->GetBinContent(i,j);
                sum4 += Histo4->GetBinContent(i,j);
                sum5 += Histo5->GetBinContent(i,j);
                sum6 += Histo6->GetBinContent(i,j);
            }
        }
        avg1 = sum1/ww/hh;
        avg2 = sum2/ww/hh;
        avg3 = sum3/ww/hh;
        avg4 = sum4/ww/hh;
        avg5 = sum5/ww/hh;
        avg6 = sum6/ww/hh;

        cout<<avg1<<'\t'<<avg2<<endl;

        for (int i = xx; i < xx+ww; i++)
        {
            for (int j = yy; j<yy+hh; j++)
            {
                Histo1->SetBinContent(i,j, avg1);
                Histo2->SetBinContent(i,j, avg2);
                Histo3->SetBinContent(i,j, avg3);
                Histo4->SetBinContent(i,j, avg4);
                Histo5->SetBinContent(i,j, avg5);
                Histo6->SetBinContent(i,j, avg6);
            }
        }
    }
// Main function
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
            for(int i = 0; i<200; i++)
            {
                for( int j=0;j<200;j++)
                {
                    DalitzNegBac->SetBinContent(i,j,DalitzNegBac->GetBinContent(i,j)/7);
                    DalitzPosBac->SetBinContent(i,j,DalitzPosBac->GetBinContent(i,j)/7);
                }
            } 
        //extract signal
            TH2F *DalitzPosSig = (TH2F*) DalitzPosCom->Clone("DalitzPosSig");
            DalitzPosSig->Add(DalitzPosBac,-1);
            TH2F *DalitzNegSig = (TH2F*) DalitzNegCom->Clone("DalitzNegSig");
            DalitzNegSig->Add(DalitzNegBac,-1);
        
        //Merge signal for exibition
        #ifdef MERGE
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0,0.25, 2.25,10);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0,0.25, 10,13);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0,0.25, 13,17);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0,0.25, 17,21);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0,0.25, 21,25.5);

            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.25,0.5, 23,26.25);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.25,0.5, 17,23);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.25,0.5, 12,17);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.25, 0.75, 5,12);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.25,0.5, 1.25,5);

            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.5,0.75, 0.75,3);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.5,0.75, 3,5);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.5,0.75, 12,19);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.5,0.75, 19,22);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.5,0.75, 22,26.5);

            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.75,1.25, 0.75,3);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.75,1.25, 3,5);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.75,1.25, 5,7);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.75,1.25, 7,10);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.75,1.25, 10,14);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.75,1.25, 14,23);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.75,1.25, 23,26.75);

            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  1.25,1.75, 1,5);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  1.25,1.75, 5,9);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  1.25,1.75, 9,15);

            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  1.75,2.75, 1.75,7);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  1.75,2.75,7,11);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  1.75,2.75, 11,15);

            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  1.25,2.25, 15,17);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  1.25,2.25, 17,21);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  1.25,2, 21,24);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  1.25,2, 24,26.5);

            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  2,3.5, 22.75,25.75);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  2,3.5, 21,22.75);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  2.25,3.5, 18.5,21);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  2.25,3.5, 15,18.5);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  2.75,4.25, 2.75,8.75);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  4.25,7.75, 4.25,7.75);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  2.75,11, 7.75,11);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  2.75,8.25,13,15);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  2.75,8.25,11,13);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  3.5,7,15,17);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  3.5,7,17,18);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  3.5,5,18,21);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  3.5,5,21,24.5);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  5,7.25,21,22.75);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  5,8,19.75,21);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  5,10,18,19.75);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  7,11.5,17,18);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  7,13,15,17);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  8.25,14,13.75,15);
            MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  8.25,13.75,11,13.75);
        #endif

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
                TCanvas *c10 = new TCanvas("c10","",CANVASIZE1,CANVASIZE2);
                DalitzPosCom->SetAxisRange(0,5,"Z");
                DalitzPosCom->Draw("colz");
                //DalitzPosCom->Draw("TEXT SAME");
                c10->SaveAs("Plots/c10_DalitzPlot_Pos_Com.pdf");

                TCanvas *c13 = new TCanvas("c13","",CANVASIZE1,CANVASIZE2);
                DalitzNegCom->SetAxisRange(0,5,"Z");
                DalitzNegCom->Draw("colz");
                //DalitzNegCom->Draw("TEXT SAME");
                c13->SaveAs("Plots/c13_DalitzPlot_Neg_Com.pdf");
            #endif

            #ifdef DRAW_BAC
                TCanvas *c14 = new TCanvas("c14","",CANVASIZE1,CANVASIZE2);
                DalitzNegBac->SetAxisRange(0,5,"Z");
                DalitzNegBac->Draw("colz");
                //DalitzNegBac->Draw("TEXT SAME");
                c14->SaveAs("Plots/c14_DalitzPlot_Neg_Bac.pdf");

                TCanvas *c11 = new TCanvas("c11","",CANVASIZE1,CANVASIZE2);
                DalitzPosBac->SetAxisRange(0,5,"Z");
                DalitzPosBac->Draw("colz");
                //DalitzPosBac->Draw("TEXT SAME");
                c11->SaveAs("Plots/c11_DalitzPlot_Pos_Bac.pdf");
            #endif


            #ifdef DRAW_SIG
                TCanvas *c12 = new TCanvas("c12","",CANVASIZE1,CANVASIZE2);
                DalitzPosSig->SetAxisRange(0,1,"Z");
                DalitzPosSig->Draw("colz");
                //DalitzPosSig->Draw("TEXT SAME");
                c12->SaveAs("Plots/c12_DalitzPlot_Pos_Sig.pdf");

                TCanvas *c15 = new TCanvas("c15","",CANVASIZE1,CANVASIZE2);
                DalitzNegSig->SetAxisRange(0,1,"Z");
                DalitzNegSig->Draw("colz");
                //DalitzNegSig->Draw("TEXT SAME");
                c15->SaveAs("Plots/c15_DalitzPlot_Pos_Sig.pdf");
            #endif

            #ifdef DRAW_ASY

                TCanvas *c16 = new TCanvas("c16","",CANVASIZE1,CANVASIZE2);
                ASSY->SetAxisRange(-0.6,0.6,"Z");
                ASSY->Draw("colz1");
                //ASSY->Draw("TEXT SAME");
                //ASSY->SetBarOffset(0.2);
                //NumMinus->SetBarOffset(-0.2);
                //NumMinus->Draw("SAME TEXT ");
                c16->SaveAs("Plots/c16_Local_Asymmytry.pdf");
            #endif
    }