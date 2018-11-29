//1 include 
    #include <TFile.h>
    #include <TH1.h>
    #include <TH2.h>
    #include <TCanvas.h>
    #include <TMath.h>
    #include <TF1.h>
    #include <iostream>
    #include "ErrorProp.hpp"
    #include "TVirtualFitter.h"
    #include <assert.h>
    using namespace std;
    #define M0B  5279.29     //Mev

//2 LIMIT
    #define LIML 5.13e3
    #define LIMH 5.7e3
    #define LIMZ 480
    #define BIN_WID 4
    #define NUM_BIN (int)((LIMH-LIML)/BIN_WID)
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
//4 choosing data file
    //#define DATA_UP
    //#define DATA_DOWN
    #define DATA_ALL

//4 Quantity of Parameters
    #define NUM_PAR_BAC 2
    #define NUM_FREE_PAR_BAC 2
    #ifdef SIGNAL_LORENTZ
        #define NUM_PAR_SIG 3
        #define NUM_FREE_PAR_SIG 3
    #endif
    #ifdef SIGNAL_CRUIJFF
        #define NUM_PAR_SIG 6
        #ifdef SIGNAL_DOUBLE_GAUS
            #define NUM_FREE_PAR_SIG 4
        #elif
            #define NUM_FREE_PAR_SIG 6
        #endif
    #endif

    #ifdef SHOULDER_GAUS
        #define NUM_PAR_SHO 3
        #define NUM_FREE_PAR_SHO 2
    #endif
    #ifdef SHOULDER_LORENTZ
        #define NUM_PAR_SHO 3
        #define NUM_FREE_PAR_SHO 2
    #endif
    #ifdef SHOULDER_ARGUS
        #define NUM_PAR_SHO 0 //TODO
    #endif
    #ifdef KAON_GAUS 
        #define NUM_PAR_KAON 3
        #define NUM_FREE_PAR_KAON 1
    #endif
    #define NUM_PAR (NUM_PAR_BAC+NUM_PAR_SIG+NUM_PAR_SHO+NUM_PAR_KAON)
    #define NUM_FREE_PAR (NUM_FREE_PAR_BAC+NUM_FREE_PAR_SHO+NUM_FREE_PAR_SIG+NUM_FREE_PAR_KAON)
    #define DEGREE_FREEDOM  (NUM_BIN - NUM_FREE_PAR)


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

//7 main function
    void fittingNeg(){
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
                TH1F *B0Neg = (TH1F*) f->Get("h_B_M0_Neg");
            #endif
            #ifdef DATA_DOWN
                TFile *f = new TFile("Output/DataMagnetDown.root");
                TH1F *B0Neg = (TH1F*) f->Get("h_B_M0_Neg");
            #endif
            #ifdef DATA_UP
                TFile *f = new TFile("Output/DataMagnetUp.root");
                TH1F *B0Neg = (TH1F*) f->Get("h_B_M0_Neg");
            #endif


        //7.1.4 Set up fitting object
        //7.1.4.1 create Neg fitting obj
            TF1 *ComNeg =      new TF1("ComNeg",   Combine,    LIML, LIMH, NUM_PAR);
    //7.2 Configure fitting parametres
        //Background
            #ifdef BACKG_EXP
                ComNeg->SetParName(NUM_PAR-NUM_PAR_KAON+1, "BackGr A");
                ComNeg->SetParameter(NUM_PAR-NUM_PAR_KAON+1, -1.13347e-03);
                ComNeg->SetParLimits(NUM_PAR-NUM_PAR_KAON+1, -1,0);

                ComNeg->SetParName(NUM_PAR-NUM_PAR_KAON+2, "BackGr B");
                ComNeg->SetParameter(NUM_PAR-NUM_PAR_KAON+2, 1.11271e+01);
                ComNeg->SetParLimits(NUM_PAR-NUM_PAR_KAON+2, 0,20);
            #endif
        //Signal
            #ifdef SIGNAL_LORENTZ 
                ComNeg->SetParName(0,"Signal I");
                ComNeg->SetParameter(0, 90);
                ComNeg->SetParLimits(0,0,400);

                ComNeg->SetParName(1,"Signal Gam");
                ComNeg->SetParameter(1, 50);
                ComNeg->SetParLimits(1,0,200);

                ComNeg->SetParName(2,"Signal x0");
                ComNeg->SetParameter(2, 5.28e3);
                ComNeg->SetParLimits(2,5200,5400);
            #endif
            #ifdef SIGNAL_CRUIJFF 
                ComNeg->SetParName(0, "Signal x0");
                ComNeg->SetParameter(0, 5.28703e+03);
                ComNeg->SetParLimits(0,5270,5290);

                ComNeg->SetParName(1, "Signal SigL");
                ComNeg->SetParameter(1, 1.91809e+01);
                ComNeg->SetParLimits(1,0,100);

                ComNeg->SetParName(2, "Signal AlpL");
                ComNeg->SetParameter(2, 0);
                ComNeg->SetParLimits(2,0,1);
                
                ComNeg->SetParName(3, "Signal SigR");
                ComNeg->SetParameter(3, 1.58689e+01);
                ComNeg->SetParLimits(3,0,100);
                
                ComNeg->SetParName(4, "Signal AlpR");
                ComNeg->SetParameter(4, 0);
                ComNeg->SetParLimits(4,0,1);
                
                ComNeg->SetParName(5, "Signal I");
                ComNeg->SetParameter(5, 1.58195e+02);
                ComNeg->SetParLimits(5,0,300);

                #ifdef SIGNAL_DOUBLE_GAUS
                ComNeg->FixParameter(2,0);
                ComNeg->FixParameter(4,0);
                #endif
            #endif
            //Shoulder
                #ifdef SHOULDER_LORENTZ
                    ComNeg->SetParName(NUM_PAR_SIG,"SHOULD I");
                    ComNeg->SetParameter(NUM_PAR_SIG, 20);
                    ComNeg->SetParLimits(NUM_PAR_SIG,0,50);

                    ComNeg->SetParName(NUM_PAR_SIG+1,"SHOULD Gam");
                    ComNeg->SetParameter(NUM_PAR_SIG+1, 100);
                    ComNeg->SetParLimits(NUM_PAR_SIG+1,5,100);

                    ComNeg->SetParName(NUM_PAR_SIG+2,"SHOULD x0");
                    ComNeg->SetParameter(NUM_PAR_SIG+2, 5100);
                    ComNeg->SetParLimits(NUM_PAR_SIG+2,4800,5100);//!THIS PARA SHOULD HAVE PHYSICS EXPLAIN
                #endif
                #ifdef SHOULDER_ARGUS
                    //TODO
                #endif
                #ifdef SHOULDER_GAUS
                    ComNeg->SetParName(NUM_PAR_SIG,"SHOULD I");
                    ComNeg->SetParameter(NUM_PAR_SIG, 5.34273e+01);
                    ComNeg->SetParLimits(NUM_PAR_SIG,0,200000);

                    ComNeg->SetParName(NUM_PAR_SIG+1,"SHOULD x0");
                    // ComNeg->SetParameter(NUM_PAR_SIG+1, 5020);
                    // ComNeg->SetParLimits(NUM_PAR_SIG+1,5020,5020);
                    ComNeg->FixParameter(NUM_PAR_SIG+1,5132);

                    ComNeg->SetParName(NUM_PAR_SIG+2,"SHOULD Sig");
                    ComNeg->SetParameter(NUM_PAR_SIG+2, 3.28215e+01);
                    ComNeg->SetParLimits(NUM_PAR_SIG+2,0,100);
                #endif
                #ifdef KAON_GAUS
                    ComNeg->SetParName(NUM_PAR_SIG+NUM_PAR_SHO,"Kaon I");
                    ComNeg->SetParameter(NUM_PAR_SIG+NUM_PAR_SHO, 0);
                    ComNeg->SetParLimits(NUM_PAR_SIG+NUM_PAR_SHO,0,0);

                    ComNeg->SetParName(NUM_PAR_SIG+NUM_PAR_SHO+1,"Kaon x0");
                    ComNeg->SetParameter(NUM_PAR_SIG+NUM_PAR_SHO+1, 5258.0107);
                    ComNeg->SetParLimits(NUM_PAR_SIG+NUM_PAR_SHO+1,5258.0107,5258.0107);

                    ComNeg->SetParName(NUM_PAR_SIG+NUM_PAR_SHO+2,"Kaon Sig");
                    ComNeg->SetParameter(NUM_PAR_SIG+NUM_PAR_SHO+2, 20);
                    ComNeg->SetParLimits(NUM_PAR_SIG+NUM_PAR_SHO+2,0,20);
                #endif
        //7.3 fitting Negtive 
            cout<<"\n*******************************POSITIVE*******************************\n"<<endl;
            TCanvas *c8 = new TCanvas("c9","",600,400);
            B0Neg->SetAxisRange(0,LIMZ,"Y");
            B0Neg->Fit(ComNeg,"R");
            double Chi2Neg = ComNeg->GetChisquare();        
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
            SigNeg->SetParErrors(&errNeg[0]);
            FourBodyNeg->SetParameters(&parNeg[NUM_PAR_SIG]);
            KaonCorrNeg->SetParameters(&parNeg[NUM_PAR_SIG+NUM_PAR_SHO]);
            BacNeg->SetParameters(&parNeg[NUM_PAR_SIG+NUM_PAR_SHO+NUM_PAR_KAON]);

            SigNeg->Draw("same");
            BacNeg->Draw("same");
            FourBodyNeg->Draw("same");
            KaonCorrNeg->Draw("SAME");

            c8->SaveAs("Plots/c9_Background&SignalFitsNeg.pdf");

    //7.5 Error Calculating
        #ifdef SIGNAL_LORENTZ
            Double_t NNNeg = parNeg[0]*parNeg[1]*TMath::Pi()/BIN_WID;
            Double_t NNNeg = parNeg[0]*parNeg[1]*TMath::Pi()/BIN_WID;
            cout<<NNNeg<<endl;
            cout<<NNNeg<<endl;
            Double_t NNeg = SigNeg->Integral(LIML,LIMH)/BIN_WID;
            Double_t NNeg = SigNeg->Integral(LIML,LIMH)/BIN_WID;
            cout<<NNeg<<endl;
            cout<<NNeg<<endl;
            Double_t ErrNNeg = ErrAMultB(NNNeg,parNeg[0], parNeg[1],errNeg[0], errNeg[1], 0);
            Double_t ErrNNeg = ErrAMultB(NNNeg,parNeg[0], parNeg[1],errNeg[0], errNeg[1],0);
        #endif

        #ifdef SIGNAL_CRUIJFF
            #ifdef SIGNAL_DOUBLE_GAUS
                Double_t NNeg = parNeg[5]*(parNeg[3]+parNeg[1])*sqrt(TMath::Pi()*2)/2/BIN_WID;
                Double_t INNeg = SigNeg->Integral(LIML,LIMH)/BIN_WID;
                Double_t IErrNNeg = ComNeg->IntegralError(M0B-55.5, M0B+55.5)/BIN_WID;
                Double_t ErrSigNeg = ErrAPlusB(errNeg[1], errNeg[3], 0);
                Double_t ErrNNeg = ErrAMultB(NNeg, parNeg[5], (parNeg[3]+parNeg[1]), errNeg[5], ErrSigNeg, 0);
            #else
            //TODO Numerical Method
            //!analytical, notgood,  Alpha tails has been ignored ,for double Gaussian it's good
            /*
                Double_t NNeg = SigNeg->Integral(LIML,LIMH);
                Double_t NNeg = SigNeg->Integral(LIML,LIMH);
                Double_t ErrNNeg = SigNeg->IntegralError(LIML,LIMH);//, parNeg, posPtr->GetCovarianceMatrix()->GetMatrixArray());
                Double_t ErrNNeg = SigNeg->IntegralError(LIML,LIMH);//, parNeg, negPtr->GetCovarianceMatrix()->GetMatrixArray());*/
                Double_t NNeg = SigNeg->Integral(LIML,LIMH)/BIN_WID;
                Double_t NNeg = SigNeg->Integral(LIML,LIMH)/BIN_WID;
                Double_t ANNeg = parNeg[5]*(parNeg[3]+parNeg[1])*sqrt(TMath::Pi()*2)/2/BIN_WID;
                Double_t ANNeg = parNeg[5]*(parNeg[3]+parNeg[1])*sqrt(TMath::Pi()*2)/2/BIN_WID;
                cout<<NNeg<<'\t'<<ANNeg<<endl;
                cout<<NNeg<<'\t'<<ANNeg<<endl;
                Double_t ErrSigNeg = ErrAPlusB(errNeg[1]/BIN_WID, errNeg[3]/BIN_WID, 0);
                Double_t ErrSigNeg = ErrAPlusB(errNeg[1]/BIN_WID, errNeg[3]/BIN_WID, 0);
                Double_t ErrNNeg = ErrAMultB(NNeg, parNeg[5]/BIN_WID, (parNeg[3]+parNeg[1])/BIN_WID, errNeg[5]/BIN_WID, ErrSigNeg, 0);
                Double_t ErrNNeg = ErrAMultB(NNeg, parNeg[5]/BIN_WID, (parNeg[3]+parNeg[1])/BIN_WID, errNeg[5]/BIN_WID, ErrSigNeg, 0);
            #endif
        #endif


    //7.6 Output results
        cout<<"\n*******************************RESULTS********************************\n"<<endl;
            cout<<"Number of B+  =\t"<<NNeg<<" +- "<<ErrNNeg<<endl;
            cout<<"INumber of B+  =\t"<<INNeg<<" +- "<<IErrNNeg<<endl;
            cout<<"Chi Squared Test+   \t"<<TMath::Prob(Chi2Neg, DEGREE_FREEDOM)<<endl;
            cout<<"Chi2 + = \t "<<Chi2Neg/DEGREE_FREEDOM<<endl;


        cout<<"\n**********************************************************************"<<endl;
        cout<<"**********************************************************************"<<endl;
        cout<<"**********************************************************************\n\n"<<endl;
    }