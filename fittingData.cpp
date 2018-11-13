#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>

#define LOWL 3800 
#define LOWR 4855
#define MIDL LOWR
#define MIDR 5110
#define HIGHL MIDR
#define HIGHR 6030

//choose Fitting method
    //#define LOWMASS_2 
    //#define LOWMASS_3
    //#define LOWMASS_4
    #define LOWMASS_EXP

    #define MIDMASS_1
    //#define MIDMASS_2

    #define HIGHMASS_3
    //#define HIGHMASS_5

    #define SIGNAL_GAUS
    //#define SIGNAL_LORENTZ

//Define BackGround fitting function
    //Lower Mass 
    double_t BackGroundLow(double_t *x, double_t *par){   
        #ifdef LOWMASS_2
        return par[0] * x[0] * x[0] + par[1] * x[0] + par[2];
        #endif
        #ifdef LOWMASS_3
        return par[0] * x[0] * x[0] * x[0] + par[1] * x[0] * x[0] + par[2] * x[0] + par[3];
        #endif
        #ifdef LOWMASS_4
        return par[0] * x[0] * x[0] * x[0] * x[0] + par[1] * x[0] * x[0] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] + par[4];
        #endif
        #ifdef LOWMASS_EXP
        return TMath::Exp(par[0] * x[0] + par[1]);
        #endif
    }
    //Med Mass BackGround
    double_t BackGroundMid(double_t *x, double_t *par)
    {
        #ifdef MIDMASS_1
        return x[0] * par[0] + par[1];
        #endif
        #ifdef MIDMASS_2
        return par[0] * x[0] * x[0] + par[1] * x[0] + par[2];
        #endif
    }
    //High mass Background
    double_t BackGroundHigh(double_t *x, double_t *par)
    {
        #ifdef HIGHMASS_3
        return par[0] * x[0] * x[0] * x[0] + par[1] * x[0] * x[0] + par[2] * x[0] + par[3];
        #endif
        #ifdef HIGHMASS_5
        return par[0] * pow(x[0],5) + par[1] * pow(x[0], 4) + par[2] * pow(x[0], 3) + par[3] * pow(x[0], 2) + par[4] * x[0] + par[5];
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
    double_t ComSigHigh(double_t *x, double_t *par)
    {
        return SignalB(x, par) + BackGroundHigh(x, &par[3]);
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
    //Lower mass fitting obj
        #ifdef LOWMASS_3
        TF1 *BGLow = new TF1("BGLow",BackGroundLow,LOWL,LOWR,4);
        BGLow->SetParameters(2e-7,7e-4,8,1e3);
        #endif
        #ifdef LOWMASS_2
        TF1 *BGLow = new TF1("BGLow",BackGroundLow,LOWL,LOWR,3);
        BGLow->SetParameters(2e-4, 1, 3e3);    
        #endif
        #ifdef LOWMASS_4
        TF1 *BGLow = new TF1("BGLow",BackGroundLow,LOWL,LOWR,5);
        BGLow->SetParameters(1e-10,6e-7,76e-3,9,1e3);
        #endif
        #ifdef LOWMASS_EXP
        TF1 *BGLow = new TF1("BGLow",BackGroundLow,LOWL,LOWR,2);
        BGLow->SetParameters(5e-3,-2);
        #endif
    //Mid mass fitting obj
        #ifdef MIDMASS_2
        TF1 *BGMid = new TF1("BGMid",BackGroundMid,MIDL,MIDR,3);
        BGMid->SetParameters(1,1,1);
        #endif
        #ifdef MIDMASS_1
        TF1 *BGMid = new TF1("BGMid",BackGroundMid,MIDL,MIDR,2);
        BGMid->SetParameters(1,1);
        #endif
    // High and Signal fitting obj
        #ifdef HIGHMASS_5
        TF1 *Com = new TF1("Com",ComSigHigh,HIGHL,HIGHR,8);
        Com->SetParameters(1,1,1,1,1,1,1,1);
        #endif
        #ifdef HIGHMASS_3
        TF1 *Com = new TF1("Com",ComSigHigh,HIGHL,HIGHR,6);
        Com->SetParameters(1,1,1,1,1,1);
        #endif
    //fitting
    B0Pos->Fit(BGLow,"R");
    B0Pos->Fit(BGMid,"R+");
    B0Pos->Fit(Com,"R+");

}





