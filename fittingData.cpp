//1 include 
    #include <TFile.h>
    #include <TH1.h>
    #include <TH2.h>
    #include <TCanvas.h>
    #include <TMath.h>
    #include <TF1.h>
    #include <iostream>
    using namespace std;

//2 LIMIT
    #define LIML 5.13e3
    #define LIMH 5.7e3
    #define LIMZ 480
    #define BIN_WID 4
    //#define LIMZ 120
//3 Choosing Interests using #define
    //3.1 choose Fitting method
        #define BACKG_EXP

        #define SIGNAL_DOUBLE_GAUS
        //#define SIGNAL_CRUIJFF
        //#define SIGNAL_LORENTZ
        
        #define SHOULDER_GAUS
        //#define SHOULDER_ARGUS
        //#define SHOULDER_LORENTZ

        #define KAON_GAUS

        //3.1.A because double Gauss and Cruijff are the same but Alpha
            #ifdef SIGNAL_DOUBLE_GAUS
            #define SIGNAL_CRUIJFF
            #endif
    //3.2 choosing data file
        //#define DATA_UP
        //#define DATA_DOWN
        #define DATA_ALL

//4 Quantity of Parameters
    #define NUM_PAR_BAC 2

    #ifdef SIGNAL_LORENTZ
        #define NUM_PAR_SIG 3
    #endif
    #ifdef SIGNAL_CRUIJFF
        #define NUM_PAR_SIG 6
    #endif

    #ifdef SHOULDER_GAUS
        #define NUM_PAR_SHO 3
    #endif
    #ifdef SHOULDER_LORENTZ
        #define NUM_PAR_SHO 3
    #endif
    #ifdef SHOULDER_ARGUS
        #define NUM_PAR_SHO 0 //TODO
    #endif
    #ifdef KAON_GAUS 
        #define NUM_PAR_KAON 3
    #endif
    #define NUM_PAR NUM_PAR_BAC+NUM_PAR_SIG+NUM_PAR_SHO+NUM_PAR_KAON


//5 Define fitting function
    //5.1 BACKGROUND_FUNCTION
        Double_t BackGround(Double_t *x, Double_t *par){   
            return TMath::Exp(par[0] * x[0] + par[1]);
        }
    //5.2 SHOULDER FUNCTION
        Double_t Shoulder(Double_t *x, Double_t *par){
            #ifdef SHOULDER_GAUS
            return par[0] * TMath::Exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));
            #endif

            #ifdef SHOULDER_LORENTZ
            return par[0]* par[1]*par[1] / TMath::Max(1.e-10, (x[0]-par[2])*(x[0]-par[2])+ par[1]*par[1]);
            #endif 

            #ifdef SHOULDER_ARGUS
            return 0 //TODO
            #endif
        }
    //5.3 Signal function
        Double_t Signal(Double_t *x, Double_t *par){

            #ifdef SIGNAL_LORENTZ
            return par[0]* par[1]*par[1] / TMath::Max(1.e-10, (x[0]-par[2])*(x[0]-par[2])+ par[1]*par[1]);
            #endif 

            #ifdef SIGNAL_CRUIJFF
                if ( x[0]<= par[0])
                    return par[5]*TMath::Exp(-(pow(x[0]-par[0],2)/(2*(pow(par[1],2)+par[2]*(pow(x[0]-par[0],2))))));
                else
                    return par[5]*TMath::Exp(-(pow(x[0]-par[0],2)/(2*(pow(par[3],2)+par[4]*(pow(x[0]-par[0],2))))));
            #endif
        }
    //Kaon correction
        Double_t KaonCorr(Double_t *x, Double_t *par)
        {
        #ifdef KAON_GAUS
            return par[0] * TMath::Exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));
        #endif
        }

    //5.4 Combine Signal and High Mass Background Function
        Double_t Combine(Double_t *x, Double_t *par){
            return Signal(x, par) + Shoulder(x,&par[NUM_PAR_SIG]) + BackGround(x, &par[NUM_PAR-NUM_PAR_KAON+1]) + KaonCorr(x,&par[NUM_PAR_SIG+NUM_PAR_SHO+NUM_PAR_KAON]);
        }
//6 Error Propagation functions
    Double_t ErrAPlusB(Double_t ErrA, Double_t ErrB)
    {
        Double_t Err2 = pow(ErrA,2)+pow(ErrB,2); 
        return TMath::Sqrt(Err2);
    }
    Double_t ErrAMinuB(Double_t ErrA, Double_t ErrB)
    {
        return ErrAPlusB(ErrA,ErrB);
    }
    Double_t ErrAMultB(Double_t C , Double_t A, Double_t B, Double_t ErrA, Double_t ErrB)
    {
        Double_t f2 = pow(C,2);
        Double_t AErrA2 = pow(ErrA/A,2);
        Double_t BErrB2 = pow(ErrB/B,2);
        Double_t Err2 = f2*(AErrA2+BErrB2);
        return TMath::Sqrt(Err2);
    }

//7 main function
    void fittingData(){
    //7.1 Declaration
        cout<<"\n\n**********************************************************************"<<endl;
        cout<<"**********************************************************************"<<endl;
        cout<<"**********************************************************************\n"<<endl;
        //7.1.1 declare fitting function
            #ifdef BACKG_EXP
                cout<<"BackGround fitting method :\tExp"<<endl;
            #endif
            #ifdef SIGNAL_CRUIJFF
                #ifdef SIGNAL_DOUBLE_GAUS
                    cout<<"Signal fitting method :\t\tDouble Gaussian"<<endl;
                #else
                    cout<<"Signal fitting method :\t\tCruijff"<<endl;
                #endif
            #endif
            #ifdef SIGNAL_LORENTZ
                cout<<"Signal fitting method :\t\tLorentz Peak"<<endl;
            #endif
            #ifdef SHOULDER_GAUS
                cout<<"Shoulder fitting method :\tGaussian"<<endl;
            #endif
            #ifdef SHOULDER_ARGUS
                cout<<"Shoulder fitting method :\tARGUS"<<endl;
            #endif
            #ifdef SHOULDER_LORENTZ
                cout<<"Shoulder fitting method :\tLorentz"<<endl;
            #endif
            #ifdef KAON_GAUS
                cout<<"Kaon Correction fitting method :\tGaussian"<<endl;
            #endif
        //7.1.2 declare Data choice
            #ifdef DATA_ALL
                cout<<"Using all Data"<<endl;
            #endif
            #ifdef DATA_UP
                cout<<"Using Magnet Up Data"<<endl;
            #endif
            #ifdef DATA_DOWN
                cout<<"Using Magnet Down Data"<<endl;
            #endif
        //7.1.3 load chosen Data from file
            #ifdef DATA_ALL
                TFile *f = new TFile("Output/DataAll.root");
                TH1F *B0Pos = (TH1F*) f->Get("h_B_M0_Pos");
                TH1F *B0Neg = (TH1F*) f->Get("h_B_M0_Neg");
            #endif
            #ifdef DATA_DOWN
                TFile *f = new TFile("Output/DataMagnetDown.root");
                TH1F *B0Pos = (TH1F*) f->Get("h_B_M0_Pos");
                TH1F *B0Neg = (TH1F*) f->Get("h_B_M0_Neg");
            #endif
            #ifdef DATA_UP
                TFile *f = new TFile("Output/DataMagnetUp.root");
                TH1F *B0Pos = (TH1F*) f->Get("h_B_M0_Pos");
                TH1F *B0Neg = (TH1F*) f->Get("h_B_M0_Neg");
            #endif


        //7.1.4 Set up fitting object
        //7.1.4.1 create Pos fitting obj
            TF1 *ComPos =      new TF1("ComPos",   Combine,    LIML, LIMH, NUM_PAR);
    //7.2 Configure fitting parametres
        //Background
            #ifdef BACKG_EXP
                ComPos->SetParName(NUM_PAR-NUM_PAR_KAON+1, "BackGr A");
                ComPos->SetParameter(NUM_PAR-NUM_PAR_KAON+1, -1.26e-3);
                ComPos->SetParLimits(NUM_PAR-NUM_PAR_KAON+1, -1,0);

                ComPos->SetParName(NUM_PAR-NUM_PAR_KAON+2, "BackGr B");
                ComPos->SetParameter(NUM_PAR-NUM_PAR_KAON+2, 10);
                ComPos->SetParLimits(NUM_PAR-NUM_PAR_KAON+2, 0,20);
            #endif
        //Signal
            #ifdef SIGNAL_LORENTZ 
                ComPos->SetParName(0,"Signal I");
                ComPos->SetParameter(0, 90);
                ComPos->SetParLimits(0,0,400);

                ComPos->SetParName(1,"Signal Gam");
                ComPos->SetParameter(1, 50);
                ComPos->SetParLimits(1,0,200);

                ComPos->SetParName(2,"Signal x0");
                ComPos->SetParameter(2, 5.28e3);
                ComPos->SetParLimits(2,5200,5400);
            #endif
            #ifdef SIGNAL_CRUIJFF 
                ComPos->SetParName(0, "Signal x0");
                ComPos->SetParameter(0, 5280);
                ComPos->SetParLimits(0,5270,5290);

                ComPos->SetParName(1, "Signal SigL");
                ComPos->SetParameter(1, 50);
                ComPos->SetParLimits(1,0,100);

                ComPos->SetParName(2, "Signal AlpL");
                ComPos->SetParameter(2, 0);
                ComPos->SetParLimits(2,0,1);
                
                ComPos->SetParName(3, "Signal SigR");
                ComPos->SetParameter(3, 50);
                ComPos->SetParLimits(3,0,100);
                
                ComPos->SetParName(4, "Signal AlpR");
                ComPos->SetParameter(4, 0);
                ComPos->SetParLimits(4,0,1);
                
                ComPos->SetParName(5, "Signal I");
                ComPos->SetParameter(5, 60);
                ComPos->SetParLimits(5,0,300);

                #ifdef SIGNAL_DOUBLE_GAUS
                ComPos->FixParameter(2,0);
                ComPos->FixParameter(4,0);
                #endif
            #endif
            //Shoulder
                #ifdef SHOULDER_LORENTZ
                    ComPos->SetParName(NUM_PAR_SIG,"SHOULD I");
                    ComPos->SetParameter(NUM_PAR_SIG, 20);
                    ComPos->SetParLimits(NUM_PAR_SIG,0,50);

                    ComPos->SetParName(NUM_PAR_SIG+1,"SHOULD Gam");
                    ComPos->SetParameter(NUM_PAR_SIG+1, 100);
                    ComPos->SetParLimits(NUM_PAR_SIG+1,5,100);

                    ComPos->SetParName(NUM_PAR_SIG+2,"SHOULD x0");
                    ComPos->SetParameter(NUM_PAR_SIG+2, 5100);
                    ComPos->SetParLimits(NUM_PAR_SIG+2,4800,5100);//!THIS PARA SHOULD HAVE PHYSICS EXPLAIN
                #endif
                #ifdef SHOULDER_ARGUS
                    //TODO
                #endif
                #ifdef SHOULDER_GAUS
                    ComPos->SetParName(NUM_PAR_SIG,"SHOULD I");
                    ComPos->SetParameter(NUM_PAR_SIG, 30);
                    ComPos->SetParLimits(NUM_PAR_SIG,0,200000);

                    ComPos->SetParName(NUM_PAR_SIG+1,"SHOULD x0");
                    ComPos->SetParameter(NUM_PAR_SIG+1, 5130);
                    ComPos->SetParLimits(NUM_PAR_SIG+1,5130,5130);

                    ComPos->SetParName(NUM_PAR_SIG+2,"SHOULD Sig");
                    ComPos->SetParameter(NUM_PAR_SIG+2, 30);
                    ComPos->SetParLimits(NUM_PAR_SIG+2,0,100);
                #endif
                #ifdef KAON_GAUS
                    ComPos->SetParName(NUM_PAR_SIG+NUM_PAR_SHO,"Kaon I");
                    ComPos->SetParameter(NUM_PAR_SIG+NUM_PAR_SHO, 0);
                    ComPos->SetParLimits(NUM_PAR_SIG+NUM_PAR_SHO,0,0);

                    ComPos->SetParName(NUM_PAR_SIG+NUM_PAR_SHO+1,"Kaon x0");
                    ComPos->SetParameter(NUM_PAR_SIG+NUM_PAR_SHO+1, 5258.0107);
                    ComPos->SetParLimits(NUM_PAR_SIG+NUM_PAR_SHO+1,5258.0107,5258.0107);

                    ComPos->SetParName(NUM_PAR_SIG+NUM_PAR_SHO+2,"Kaon Sig");
                    ComPos->SetParameter(NUM_PAR_SIG+NUM_PAR_SHO+2, 10);
                    ComPos->SetParLimits(NUM_PAR_SIG+NUM_PAR_SHO+2,10,10);
                #endif
        //7.4.3 Create Neg fitting obj
            TF1 *ComNeg = ComPos;
            ComNeg->SetName("ComNeg");
    //7.3 fitting Postive 
        cout<<"\n*******************************POSITIVE*******************************\n"<<endl;
        TCanvas *c8 = new TCanvas("c8","",600,400);
        B0Pos->SetAxisRange(0,LIMZ,"Y");
        B0Pos->Fit(ComPos,"R");
        double Chi2Pos = ComPos->GetChisquare();
        //B0Pos->SetAxisRange(5100,5500);
        
        TF1 *BacPos      = new TF1("Bac",   BackGround, LIML, LIMH, NUM_PAR_BAC);
        TF1 *FourBodyPos = new TF1("4Body", Shoulder,   LIML, LIMH, NUM_PAR_SHO); 
        TF1 *SigPos      = new TF1("Sig",   Signal,     LIML, LIMH, NUM_PAR_SIG);
        TF1 *KaonCorrPos = new TF1("KaonCorrPos",KaonCorr,LIML, LIMH, NUM_PAR_KAON);
        SigPos->SetLineColor(kBlue);
        BacPos->SetLineColor(kYellow);
        FourBodyPos->SetLineColor(kGreen);
        KaonCorrPos->SetLineColor(kMagenta);

        Double_t parPos[NUM_PAR];
        Double_t *errPos;
        ComPos->GetParameters(parPos);
        errPos = (Double_t*)ComPos->GetParErrors();
        SigPos->SetParameters(&parPos[0]);
        FourBodyPos->SetParameters(&parPos[NUM_PAR_SIG]);
        KaonCorrPos->SetParameters(&parPos[NUM_PAR_SIG+NUM_PAR_SHO]);
        BacPos->SetParameters(&parPos[NUM_PAR - NUM_PAR_BAC+1]);

        SigPos->Draw("same");
        BacPos->Draw("same");
        FourBodyPos->Draw("same");
        KaonCorrPos->Draw("SAME");

        c8->SaveAs("Plots/c8_Background&SignalFitsPos.pdf");

    //7.4 fitting Negtive
        cout<<"\n*******************************NEGTIVE*******************************\n"<<endl;
        TCanvas *c9 = new TCanvas("c9","",600,400);
        B0Neg->SetAxisRange(0,LIMZ,"Y");
        B0Neg->Fit(ComNeg,"R");
        double Chi2Neg = ComNeg->GetChisquare();
        //B0Neg->SetAxisRange(5100,5500);

        TF1 *BacNeg      = new TF1("Bac",   BackGround, LIML, LIMH, NUM_PAR_BAC);
        TF1 *FourBodyNeg = new TF1("4Body", Shoulder,   LIML, LIMH, NUM_PAR_SHO); 
        TF1 *SigNeg      = new TF1("Sig",   Signal,     LIML, LIMH, NUM_PAR_SIG);
        TF1 *KaonCorrNeg = new TF1("KaonCorrNeg",KaonCorr,LIML, LIMH, NUM_PAR_KAON);
        SigNeg->SetLineColor(kBlue);
        BacNeg->SetLineColor(kYellow);
        FourBodyNeg->SetLineColor(kGreen);
        KaonCorrNeg->SetLineColor(kMagenta);

        Double_t parNeg[NUM_PAR];
        Double_t *errNeg;
        ComNeg->GetParameters(parNeg);
        errNeg = (Double_t*)ComNeg->GetParErrors();
        SigNeg->SetParameters(&parNeg[0]);
        FourBodyNeg->SetParameters(&parNeg[NUM_PAR_SIG]);
        KaonCorrNeg->SetParameters(&parPos[NUM_PAR_SIG+NUM_PAR_SHO]);
        BacNeg->SetParameters(&parNeg[NUM_PAR-NUM_PAR_BAC+1]);
        

        SigNeg->Draw("same");
        BacNeg->Draw("same");
        FourBodyNeg->Draw("same");
        KaonCorrPos->Draw("SAME");

        c9->SaveAs("Plots/c9_Background&SignalFitsNeg.pdf");

    //7.5 Error Calculating
        #ifdef SIGNAL_LORENTZ
            Double_t NNPos = parPos[0]*parPos[1]*TMath::Pi()/BIN_WID;
            Double_t NNNeg = parNeg[0]*parNeg[1]*TMath::Pi()/BIN_WID;
            cout<<NNPos<<endl;
            cout<<NNNeg<<endl;
            Double_t NPos = SigPos->Integral(LIML,LIMH)/BIN_WID;
            Double_t NNeg = SigNeg->Integral(LIML,LIMH)/BIN_WID;
            cout<<NPos<<endl;
            cout<<NNeg<<endl;
            Double_t ErrNPos = ErrAMultB(NNPos,parPos[0]/BIN_WID, parPos[1]/BIN_WID,errPos[0]/BIN_WID, errPos[1]/BIN_WID);
            Double_t ErrNNeg = ErrAMultB(NNNeg,parNeg[0]/BIN_WID, parNeg[1]/BIN_WID,errNeg[0]/BIN_WID, errNeg[1]/BIN_WID);
        #endif

        #ifdef SIGNAL_CRUIJFF
            #ifdef SIGNAL_DOUBLE_GAUS
                Double_t NPos = SigPos->Integral(LIML,LIMH)/BIN_WID;
                Double_t NNeg = SigNeg->Integral(LIML,LIMH)/BIN_WID;
                Double_t ANPos = parPos[5]*(parPos[3]+parPos[1])*sqrt(TMath::Pi()*2)/2/BIN_WID;
                Double_t ANNeg = parNeg[5]*(parNeg[3]+parNeg[1])*sqrt(TMath::Pi()*2)/2/BIN_WID;
                cout<<NPos<<'\t'<<ANPos<<endl;
                cout<<NNeg<<'\t'<<ANNeg<<endl;
                Double_t ErrSigPos = ErrAPlusB(errPos[1]/BIN_WID, errPos[3]/BIN_WID);
                Double_t ErrSigNeg = ErrAPlusB(errNeg[1]/BIN_WID, errNeg[3]/BIN_WID);
                Double_t ErrNPos = ErrAMultB(NPos, parPos[5]/BIN_WID, (parPos[3]+parPos[1])/BIN_WID, errPos[5]/BIN_WID, ErrSigPos);
                Double_t ErrNNeg = ErrAMultB(NNeg, parNeg[5]/BIN_WID, (parNeg[3]+parNeg[1])/BIN_WID, errNeg[5]/BIN_WID, ErrSigNeg);
            #else
            //TODO Numerical Method
            //!analytical, notgood,  Alpha tails has been ignored ,for double Gaussian it's good
            /*
                Double_t NPos = SigPos->Integral(LIML,LIMH);
                Double_t NNeg = SigNeg->Integral(LIML,LIMH);
                Double_t ErrNPos = SigPos->IntegralError(LIML,LIMH);//, parPos, posPtr->GetCovarianceMatrix()->GetMatrixArray());
                Double_t ErrNNeg = SigNeg->IntegralError(LIML,LIMH);//, parNeg, negPtr->GetCovarianceMatrix()->GetMatrixArray());*/
                Double_t NPos = SigPos->Integral(LIML,LIMH)/BIN_WID;
                Double_t NNeg = SigNeg->Integral(LIML,LIMH)/BIN_WID;
                Double_t ANPos = parPos[5]*(parPos[3]+parPos[1])*sqrt(TMath::Pi()*2)/2/BIN_WID;
                Double_t ANNeg = parNeg[5]*(parNeg[3]+parNeg[1])*sqrt(TMath::Pi()*2)/2/BIN_WID;
                cout<<NPos<<'\t'<<ANPos<<endl;
                cout<<NNeg<<'\t'<<ANNeg<<endl;
                Double_t ErrSigPos = ErrAPlusB(errPos[1],errPos[3]);
                Double_t ErrSigNeg = ErrAPlusB(errNeg[1],errNeg[3]);
                Double_t ErrNPos = ErrAMultB(NPos, parPos[5], parPos[3]+parPos[1], errPos[5], ErrSigPos);
                Double_t ErrNNeg = ErrAMultB(NNeg, parNeg[5], parNeg[3]+parNeg[1], errNeg[5], ErrSigNeg);
            #endif
        #endif
        Double_t Asym = (NNeg-NPos)/(NNeg+NPos);
        Double_t ErrN    = ErrAPlusB(ErrNPos, ErrNNeg);
        Double_t EffPos  = NPos/(NPos+NNeg);
        Double_t EffNeg  = NNeg/(NPos+NNeg);

        Double_t ErrEffPos = ErrAMultB(EffPos,NPos, NPos+NNeg, ErrNPos,ErrN);
        Double_t ErrEffNeg = ErrAMultB(EffNeg,NNeg, NPos+NNeg, ErrNNeg,ErrN);

        Double_t ErrA      = ErrAMinuB(ErrEffNeg, ErrEffPos);

    //7.6 Output results
        cout<<"\n*******************************RESULTS********************************\n"<<endl;
        cout<<"Number of B+  =\t"<<NPos<<" +- "<<ErrNPos<<endl;
        cout<<"Number of B-  =\t"<<NNeg<<" +- "<<ErrNNeg<<endl;
        cout<<"Efficiency N+ =\t"<<EffPos<<" +- "<<ErrEffPos<<endl;
        cout<<"Efficiency N- =\t"<<EffNeg<<" +- "<<ErrEffNeg<<endl;
        cout<<"Effi Glb Asym =\t"<<Asym<<" +- "<<ErrA<<endl;
        ErrA = ErrAMultB(Asym, NNeg-NPos, NNeg+NPos, ErrN, ErrN);
        cout<<"Devi Glb Asym =\t"<<Asym<<" +- "<<ErrA<<endl;
        cout<<"Predicted Error =\t"<<TMath::Sqrt((1-Asym*Asym)/(NPos+NNeg))<<endl;
        cout<<"Chi2 Pos and Neg\t "<<Chi2Pos/8<<'\t'<<Chi2Neg/8<<endl;
        cout<<"\n**********************************************************************"<<endl;
        cout<<"**********************************************************************"<<endl;
        cout<<"**********************************************************************\n\n"<<endl;
    }