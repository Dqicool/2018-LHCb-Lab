#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <iostream>
using namespace std;

#define LIML 5.05e3
#define LIMH 6e3

//choose Fitting method
    //#define MASS_2 
    //#define MASS_3
    //#define MASS_4
    #define MASS_EXP


    //#define SIGNAL_GAUS
    #define SIGNAL_LORENTZ

//Define BackGround fitting function
    //Lower Mass 
    Double_t BackGround(Double_t *x, Double_t *par){   
        #ifdef MASS_2
        return par[0] * pow(x[0] - par[1], 2) - par[2];
        #endif
        #ifdef MASS_3
        return par[0] * x[0] * x[0] * x[0] + par[1] * x[0] * x[0] + par[2] * x[0] + par[3];
        #endif
        #ifdef MASS_4
        return par[0] * x[0] * x[0] * x[0] * x[0] + par[1] * x[0] * x[0] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] + par[4];
        #endif
        #ifdef MASS_EXP
        return TMath::Exp(par[0] * x[0] + par[1]);
        #endif
    }
    
//Define Signal fitting function
    Double_t Signal(Double_t *x, Double_t *par)
    {
        #ifdef SIGNAL_GAUS
        return TMath::Max(par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2)),1e-4);
        #endif
        #ifdef SIGNAL_LORENTZ
        return par[0]* par[1]*par[1] / TMath::Max(1.e-10, (x[0]-par[2])*(x[0]-par[2])+ par[1]*par[1]);
        #endif
    }
//Combine Signal and High Mass Background Function
     Double_t Combine(Double_t *x, Double_t *par)
    {
        return Signal(x, par) + Signal(x,&par[3]) + BackGround(x, &par[6]);
    }
//Error Propagation functions
    Double_t ErrAPlusB(Double_t ErrA, Double_t ErrB)
    {
        Double_t Err2 = pow(ErrA,2)+pow(ErrB,2); 
        return sqrt(Err2);
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
        return sqrt(Err2);
    }
//main function
void fittingData(){
//Loading files
    TFile *f1 = new TFile("Output/DataAll.root");
    TFile *f2 = new TFile("Output/DataMagnetDown.root");
    TFile *f3 = new TFile("Output/DataMagnetUp.root");
    TFile *f4 = new TFile("Output/PhaseSpace.root");
//Define needed Data from file
    TH1F *B0Pos = (TH1F*) f1->Get("h_B_M0_Pos");
    TH1F *B0Neg = (TH1F*) f1->Get("h_B_M0_Neg");
//Set up function object
    //fitting obj
        #ifdef MASS_2
        TF1 *Com = new TF1("Com",Combine,LIML,LIMH,6);
        Com->SetParameters(90,50,5.280e3,3.9e-4,5440,20);    
        TF1 *Bac = new TF1("Bac",BackGround,LIML,LIMH,3);
        TF1 *Sig = new TF1("Sig",Signal,LIML,LIMH,3);
        #endif
        
        #ifdef MASS_3
        TF1 *Com = new TF1("Com",Combine,LIML,LIMH,7);
        Com->SetParameters(90,50,5.280e3,1e3);
        TF1 *Bac = new TF1("Bac",BackGround,LIML,LIMH,4);
        TF1 *Sig = new TF1("Sig",Signal,LIML,LIMH,3);
        #endif
        
        #ifdef MASS_4
        TF1 *Com = new TF1("Com",Combine,LIML,LIMH,8);
        Com->SetParameters(90,50,5.280e3,1e-9,1e-6,1e-4,1e-1,1e3);
        TF1 *Bac = new TF1("Bac",BackGround,LIML,LIMH,5);
        TF1 *Sig = new TF1("Sig",Signal,LIML,LIMH,3);
        #endif
        
        #ifdef MASS_EXP
        TF1 *Com = new TF1("Com",Combine,LIML,LIMH,8);
        Com->SetParameters(90,50,5.280e3,  20,100,5100,  -0.0025,17);
	
        TF1 *Bac = new TF1("Bac",BackGround,LIML,LIMH,2);
        TF1 *Sig = new TF1("Sig",Signal,LIML,LIMH,3);
	    TF1 *FourBody = new TF1("4Body",Signal,LIML,LIMH,3); 
        #endif

    //fitting
        //Setting Parametar limits
        Com->SetParLimits(0,0,100);
        Com->SetParLimits(1,0,200);
        Com->SetParLimits(2,5200,5400);
        Com->SetParLimits(3,0,50);
        Com->SetParLimits(4,5,100);
        Com->SetParLimits(5,4800,5100);
        Com->SetParLimits(7,0,1e5);

        //Naming Parameters
        Com->SetParNames(   "Signal I", "Signal Gam", "Signal x0",
                            "4 Body I", "4 Body Gam", "4 body x0",
                            "  BacG A",   "  BacG B");
    
        //Postive 
            TCanvas *c8 = new TCanvas("c8","",600,400);
            B0Pos->SetAxisRange(0,120,"Y");
            B0Pos->Fit(Com,"R");

            Double_t ErrPosI = Com->GetParError(0);
            Double_t ErrPosGamma = Com->GetParError(1);
            
            Double_t parPos[20];
            Com->GetParameters(parPos);
            Sig->SetParameters(&parPos[0]);
            FourBody->SetParameters(&parPos[3]);
            Bac->SetParameters(&parPos[6]);
            
            Sig->SetLineColor(kBlue);
            Bac->SetLineColor(kYellow);
            FourBody->SetLineColor(kGreen);

            Sig->Draw("same");
            Bac->Draw("same");
            FourBody->Draw("same");
            c8->SaveAs("Plots/c8_Background&SignalFitsPos.pdf");

        //Negtive
            TCanvas *c9 = new TCanvas("c9","",600,400);
            B0Neg->SetAxisRange(0,120,"Y");
            B0Neg->Fit(Com,"R");

            Double_t ErrNegI = Com->GetParError(0);
            Double_t ErrNegGamma = Com->GetParError(1);

            Double_t parNeg[20];
            Com->GetParameters(parNeg);
            Sig->SetParameters(&parNeg[0]);
            FourBody->SetParameters(&parNeg[3]);
            Bac->SetParameters(&parNeg[6]);
            
            Sig->SetLineColor(kBlue);
            Bac->SetLineColor(kYellow);
            FourBody->SetLineColor(kGreen);
            Sig->Draw("same");
            Bac->Draw("same");
            FourBody->Draw("same");
            c9->SaveAs("Plots/c9_Background&SignalFitsNeg.pdf");

    //Error Calculating
        Double_t NPos = parPos[0]*parPos[1]*TMath::Pi();
        Double_t NNeg = parNeg[0]*parNeg[1]*TMath::Pi();
        Double_t Asym = (NNeg-NPos)/(NNeg+NPos);
        Double_t ErrNPos = ErrAMultB(NPos,parPos[0], parPos[1],ErrPosI, ErrPosGamma);
        Double_t ErrNNeg = ErrAMultB(NNeg,parNeg[0], parNeg[1],ErrNegI, ErrNegGamma);
        Double_t ErrN    = ErrAPlusB(ErrNPos, ErrNNeg);
        Double_t EffPos  = NPos/(NPos+NNeg);
        Double_t EffNeg  = NNeg/(NPos+NNeg);
        Double_t ErrEffPos = ErrAMultB(EffPos,NPos, NPos+NNeg, ErrNPos,ErrN);
        Double_t ErrEffNeg = ErrAMultB(EffNeg,NNeg, NPos+NNeg, ErrNNeg,ErrN);
        Double_t ErrA      = ErrAMinuB(ErrEffNeg, ErrEffPos);

        ErrA = ErrAMultB(Asym, NNeg-NPos, NNeg+NPos, ErrN, ErrN);
    //Output

        cout<<"Number of B+  =\t"<<NPos<<" +- "<<ErrNPos<<endl;
        cout<<"Number of B-  =\t"<<NNeg<<" +- "<<ErrNNeg<<endl;
        cout<<"Efficiency N+ =\t"<<EffPos<<" +- "<<ErrEffPos<<endl;
        cout<<"Efficiency N- =\t"<<EffNeg<<" +- "<<ErrEffNeg<<endl;
        cout<<"Global Asymme =\t"<<Asym<<" +- "<<ErrA<<endl;
        cout<<"Predicted Error =\t"<<sqrt((1-Asym*Asym)/(NPos+NNeg))<<endl;
}



