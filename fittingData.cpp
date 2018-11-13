#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>

#define LIML 5.11e3
#define LIMH 5.5e3


//choose Fitting method
    //#define MASS_2 
    //#define MASS_3
    //#define MASS_4
    #define MASS_EXP


    #define SIGNAL_GAUS
    //#define SIGNAL_LORENTZ

//Define BackGround fitting function
    //Lower Mass 
    double_t BackGround(double_t *x, double_t *par){   
        #ifdef MASS_2
        return par[0] * x[0] * x[0] + par[1] * x[0] + par[2];
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
    double_t SignalB(double_t *x, double_t *par)
    {
        #ifdef SIGNAL_GAUS
        return TMath::Max(par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2)),1e-4);
        #endif
        #ifdef SIGNAL_LORENTZ
        return TMath::Max(1e-4,(0.5*par[0]*par[1]/TMath::Pi()) / TMath::Max(1.e-10, (x[0]-par[2])*(x[0]-par[2])+ 0.25*par[1]*par[1]));
        #endif
    }
//Combine Signal and High Mass Background Function
    double_t Combine(double_t *x, double_t *par)
    {
        return SignalB(x, par) + BackGround(x, &par[3]);
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
    TH1F *B0Neg = (TH1F*) f2->Get("h_B_M0_Neg");
//Set up function object
    //fitting obj
        #ifdef MASS_3
        TF1 *Com = new TF1("Com",Combine,LIML,LIMH,6);
        Com->SetParameters(2e-7,7e-4,8,1e3);
        #endif
        #ifdef MASS_2
        TF1 *Com = new TF1("Com",Combine,LIML,LIMH,5);
        Com->SetParameters(2e-4, 1, 3e3);    
        #endif
        #ifdef MASS_4
        TF1 *Com = new TF1("Com",Combine,LIML,LIMH,7);
        Com->SetParameters(1e-10,6e-7,76e-3,9,1e3);
        #endif
        #ifdef MASS_EXP
        TF1 *Com = new TF1("Com",Combine,LIML,LIMH,5);
        Com->SetParameters(5e-3,-2);
        #endif

    //fitting

    B0Pos->Fit(Com,"R+");

}





