//1 include 
    #include <TFile.h>
    #include <TH1.h>
    #include <TH2.h>
    #include <TCanvas.h>
    #include <TMath.h>
    #include <TF1.h>
    #include <TF2.h>
    #include <iostream>
//Canvasize
    using namespace std;
    #define CANVASIZE1 600
    #define CANVASIZE2 400
// Draw select
    #define DRAW_BAC
    #define DRAW_COM
    #define DRAW_SIG
    #define DRAW_ASY
    #define DRAW_SIGNIF
// Merge select
    #define MANUAL_MERGE
    //#define AUTO_MERGE
// Define Data Choice
    #define DATA_ALL
    //#define DATA_UP
    //#define DATA_DOWN
    //#define DATA_PHASE_SPACE
// Merge func 
    void MergeBins(TH2F* Histo1, TH2F* Histo2, TH2F* Histo3, TH2F* Histo4,TH2F* Histo5,TH2F* Histo6, double x, double w, double y, double h)
    {
        int xx = 4*x+1; int yy = 4*y+1; int ww = 4*(w-x); int hh = 4*(h-y);
        Double_t sum1 = 0, avg1;
        Double_t sum2 = 0, avg2;
        Double_t sum3 = 0, avg3;
        Double_t sum4 = 0, avg4;
        Double_t sum5 = 0, avg5;
        Double_t sum6 = 0, avg6;
        Double_t sumErr1 = 0, avgErr1;
        Double_t sumErr2 = 0, avgErr2;
        Double_t sumErr3 = 0, avgErr3;
        Double_t sumErr4 = 0, avgErr4;
        Double_t sumErr5 = 0, avgErr5;
        Double_t sumErr6 = 0, avgErr6;
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
                sumErr1+= pow(Histo1->GetBinError(i,j),2);
                sumErr2+= pow(Histo2->GetBinError(i,j),2);
                sumErr3+= pow(Histo3->GetBinError(i,j),2);
                sumErr4+= pow(Histo4->GetBinError(i,j),2);
                sumErr5+= pow(Histo5->GetBinError(i,j),2);
                sumErr6+= pow(Histo6->GetBinError(i,j),2);
            }
        }
        avg1 = sum1/ww/hh;
        avg2 = sum2/ww/hh;
        avg3 = sum3/ww/hh;
        avg4 = sum4/ww/hh;
        avg5 = sum5/ww/hh;
        avg6 = sum6/ww/hh;
        avgErr1 = sqrt(sumErr1)/ww/hh;
        avgErr2 = sqrt(sumErr2)/ww/hh;
        avgErr3 = sqrt(sumErr3)/ww/hh;
        avgErr4 = sqrt(sumErr4)/ww/hh;
        avgErr5 = sqrt(sumErr5)/ww/hh;
        avgErr6 = sqrt(sumErr6)/ww/hh;

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

                Histo1->SetBinError(i,j, avgErr1);
                Histo2->SetBinError(i,j, avgErr2);
                Histo3->SetBinError(i,j, avgErr3);
                Histo4->SetBinError(i,j, avgErr4);
                Histo5->SetBinError(i,j, avgErr5);
                Histo6->SetBinError(i,j, avgErr6);
            }
        }
    }
     void MergeBins1(TH2F* Histo1, double x, double w, double y, double h)
    {
        int xx = x+1; int yy = y+1; int ww = (w-x); int hh = (h-y);
        Double_t sum1 = 0, avg1;
        Double_t sumErr1 = 0, avgErr1;
        for (int i = xx; i < xx+ww; i++)
        {
            for (int j = yy; j<yy+hh; j++)
            {
                sum1 += Histo1->GetBinContent(i,j);
                sumErr1+= pow(Histo1->GetBinError(i,j),2);
            }
        }
        avg1 = sum1/ww/hh;
        avgErr1 = sqrt(sumErr1)/ww/hh;

        for (int i = xx; i < xx+ww; i++)
        {
            for (int j = yy; j<yy+hh; j++)
            {
                Histo1->SetBinContent(i,j, avg1);

                Histo1->SetBinError(i,j, avgErr1);
            }
        }
    }
//Error Getting function
    TH2F* GetHistErr(TH2F* Histo)
    {
        TH2F* ErrHist = (TH2F*)Histo->Clone();
        for(int i = 0;i<200;i++)
        {
            for(int j = 0; j<200;j++)
            {
                ErrHist->SetBinContent(i,j,Histo->GetBinError(i,j));
            }
        }
        return ErrHist;
    }
//absolute value getting func
    TH2F* GetHistAbs(TH2F* Histo)
    {
        TH2F* AbsHist = (TH2F*)Histo->Clone();
        for(int i = 0;i<200;i++)
        {
            for(int j = 0; j<200;j++)
            {
                double tmp;
                if (Histo->GetBinContent(i,j)<0)
                    tmp= -Histo->GetBinContent(i,j);
                else 
                    tmp = Histo->GetBinContent(i,j);
                AbsHist->SetBinContent(i,j,tmp);
            }
        }
        return AbsHist;
    }
//histo Histo partial sum function
    double SumPartial(TH2F* Histo, int xll, int yll, int xul, int yul)
    {
        double res = 0;
        for(int i = xll; i<=xul;i++)
        {
            for(int j = yll; j<=yul;j++)
            {
                res+= Histo->GetBinContent(i,j);
            }
        }
        return res;
    }


//Auto Merger function
    void AutoMerge(TH2F *Signal, int xll, int yll, int xul, int yul, double NumPerCell)
    {
        if (yll<xll)
            yll = xll;
        if (yul>112-xll)
            yul = 112-xll;
        if (xul>56)
            xul=56;
        double NumAll   = SumPartial(Signal, xll,        yll,        xul,      yul);
        double DiffUD, DiffLR, DevideLR, DevideUD;
        if (xul - xll<2 && yul - yll<2)
        {
            return;
        }
        else if (xul - xll<2)
        {
            DiffLR = 0;
            DevideUD = yll+(yul-yll)/2;

            double NumDown  = SumPartial(Signal, xll,        yll,        xul,      DevideUD);
            double NumUp    = SumPartial(Signal, xll,        DevideUD+1, xul,      yul);

            DiffUD = TMath::Abs(NumUp - NumDown);
        }
        else if (yul - yll<2)
        {
            DiffUD = 0;
            DevideLR = xll+(xul-xll)/2;

            double NumLeft  = SumPartial(Signal, xll,        yll,        DevideLR, yul);
            double NumRight = SumPartial(Signal, DevideLR+1, yll,        xul,      yul);

            DiffLR = TMath::Abs(NumLeft - NumRight);
        }
        else
        {
            DevideLR = xll+(xul-xll)/2;
            DevideUD = yll+(yul-yll)/2;

            double NumLeft  = SumPartial(Signal, xll,        yll,        DevideLR, yul);
            double NumRight = SumPartial(Signal, DevideLR+1, yll,        xul,      yul);
            double NumDown  = SumPartial(Signal, xll,        yll,        xul,      DevideUD);
            double NumUp    = SumPartial(Signal, xll,        DevideUD+1, xul,      yul);

            DiffUD = TMath::Abs(NumUp - NumDown);
            DiffLR = TMath::Abs(NumLeft - NumRight);
        }
        if (NumAll<=NumPerCell)
        {
            cout<<"UD\t"<<NumAll<<endl;
            MergeBins1(Signal,xll, xul, yll,yul);
            //cout<<"MERGE\t"<<(double)xll/4<<'\t'<<(double)xul/4<<'\t'<<(double)yll/4<<'\t'<<(double)yul/4<<'\t'<<NumAll<<endl;
            //cout<<"MIN END"<<endl;
            return;
        }

        if(DiffLR>DiffUD)
        {
            //cout<<"LR\t"<<(double)xll/4<<'\t'<<(double)xul/4<<'\t'<<(double)yll/4<<'\t'<<(double)yul/4<<'\t'<<NumAll<<'\t'<<DiffLR<<'\t'<<DiffUD<<endl;
            
            AutoMerge(Signal, xll,        yll,        DevideLR, yul, NumPerCell);
            AutoMerge(Signal, DevideLR, yll,        xul,      yul, NumPerCell);
            return;
        }
        else
        {
            //cout<<"ud\t"<<(double)xll/4<<'\t'<<(double)xul/4<<'\t'<<(double)yll/4<<'\t'<<(double)yul/4<<'\t'<<NumAll<<'\t'<<(double)DevideUD/4<<endl;
            AutoMerge(Signal, xll,        yll,        xul,      DevideUD, NumPerCell);
            AutoMerge(Signal, xll,        DevideUD,   xul,      yul,      NumPerCell);
            return;
        }
    }
//Chi2 test function
    double BinnedChi2(TH2F *Pos, TH2F *Neg, int &ndf, int xll, int xhl, int yll, int yhl)
    {
        TH2F * Scp = (TH2F*)Neg->Clone("Scp");
        double alpha = SumPartial(Pos, xll,yll,xhl,yhl)/SumPartial(Neg, xll,yll,xhl,yhl);
        //cout<<alpha<<endl;
        ndf = -1;
        for (int i = xll;i<xhl;i++)
        {
            for (int j=yll;j<yhl;j++)
            {
                if (Pos->GetBinContent(i,j)== 0 || Neg->GetBinContent(i,j)== 0){
                    Scp->SetBinContent(i,j,0);
                }
                else{
                    //cout<<Pos->GetBinContent(i,j)<<'\t'<<Neg->GetBinContent(i,j)<<endl;
                    double tmp = (Pos->GetBinContent(i,j)-alpha*Neg->GetBinContent(i,j))/sqrt(Pos->GetBinContent(i,j)+alpha*alpha*Neg->GetBinContent(i,j));
                    cout<<Pos->GetBinContent(i,j)<<'\t'<<Neg->GetBinContent(i,j)<<'\t'<<tmp<<endl;
                    Scp->SetBinContent(i,j,tmp);
                    ndf++;
                }
            }
        }
        double chi2 = 0;
        for (int i = xll;i<xhl;i++)
        {
            for (int j=yll;j<yhl;j++)
            {
                chi2 += pow(Scp->GetBinContent(i,j),2);
            }
        }
        return chi2;
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
        
        //Merging signal for exibition manually
            #ifdef MANUAL_MERGE
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0,0.25, 2.25,10);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0,0.25, 10,13);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0,0.25, 13,17);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0,0.25, 17,21);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0,0.25, 21,25.5);

                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.25,0.5, 23,26.5);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.25,0.5, 17,23);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.25,0.5, 12,17);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.25, 0.75, 5,12);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  0.25,0.5, 1,5);

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
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  3.5,8.25,15,17);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  3.5,7,17,18);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  3.5,5,18,21);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  3.5,5,21,24.5);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  5,7.25,21,23);
                //MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  5,8,20,21);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  5,10,18,21);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  7,11.5,17,18);

                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  7,13,15,17);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  8.25,14,13.75,15);
                MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  8.25,13.75,11,13.75);
                //MergeBins(DalitzPosSig, DalitzNegSig, DalitzPosBac, DalitzPosCom, DalitzNegBac, DalitzNegCom,  8.25,14,11,17);

            #endif

        //7.5 Calculate Global asymmytries
            TH2F *NumAll = (TH2F*) DalitzPosSig->Clone("NumAll");
            NumAll->Add(DalitzNegSig);
            Double_t NumB = 0;
            for(int i = 0; i<200; i++)
            {
                for( int j=0;j<200;j++)
                {
                    NumB += NumAll->GetBinContent(i,j);
                }
            } 

            TH2F *NumMinus = (TH2F*)DalitzNegSig->Clone("NumMinus");
            NumMinus->Add(DalitzPosSig,-1);
            Double_t NumA=0;
            for(int i = 0; i<200; i++)
            {
                for( int j=0;j<200;j++)
                {
                    NumA += NumMinus->GetBinContent(i,j);
                }
            } 
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
                DalitzPosSig->SetAxisRange(0,5,"Z");
                DalitzPosSig->Draw("colz");
                //DalitzPosSig->Draw("TEXT SAME");
                c12->SaveAs("Plots/c12_DalitzPlot_Pos_Sig.pdf");

                TCanvas *c15 = new TCanvas("c15","",CANVASIZE1,CANVASIZE2);
                DalitzNegSig->SetAxisRange(0,5,"Z");
                DalitzNegSig->Draw("colz");
                //DalitzNegSig->Draw("TEXT SAME");
                c15->SaveAs("Plots/c15_DalitzPlot_Pos_Sig.pdf");
            #endif

            #ifdef DRAW_ASY

                TCanvas *c16 = new TCanvas("c16","",CANVASIZE1,CANVASIZE2);
                ASSY->SetAxisRange(-1,1,"Z");
                ASSY->Draw("colz1");
                //ASSY->Draw("TEXT SAME");
                //ASSY->SetBarOffset(0.2);
                //NumMinus->SetBarOffset(-0.2);
                //NumMinus->Draw("SAME TEXT ");
                c16->SaveAs("Plots/c16_Local_Asymmytry.pdf");
            #endif
            #ifdef DRAW_SIGNIF
                TCanvas *c17 = new TCanvas("c17","",CANVASIZE1,CANVASIZE2);
                TH2F* Err = GetHistErr(ASSY);
                TH2F* SIGNIF = (TH2F*)ASSY->Clone("absSignif of Asymm");
                SIGNIF->SetAxisRange(0,5,"Z");
                TH2F *absSIGNIF = GetHistAbs(SIGNIF);
                absSIGNIF->Divide(Err);
                absSIGNIF->Draw("colz");
                Err->SetAxisRange(0,1,"Z");
                //Err->Draw("colz1");
                c17->SaveAs("Plots/c17_Local_Asymmytry_absSignif.pdf");

                TCanvas *c18 = new TCanvas("c18","",CANVASIZE1,CANVASIZE2);
                SIGNIF->Divide(Err);
                SIGNIF->Draw("colz1");
                SIGNIF->SetAxisRange(-5,5,"Z");
                c18->SaveAs("Plots/c18_Local_Asymmytry_Signif.pdf");
            #endif
            #ifdef AUTO_MERGE
                #define NUM_PER_CELL 100
                    TCanvas *c18 = new TCanvas("c18","",CANVASIZE1,CANVASIZE2);
                    AutoMerge(DalitzNegSig,0,0,60,120,NUM_PER_CELL);
                    //DalitzNegSig->SetAxisRange(0,20,"Z");
                    DalitzNegSig->Draw("colz");
            #endif
            int ndf;
            double chi2 = BinnedChi2(DalitzPosSig,DalitzNegSig,ndf,1*4,3*4,5*4, 15*4);
            cout<<chi2<<'\t'<<ndf<<endl;
            cout<<TMath::Prob(chi2,ndf)<<endl;
    }