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

    double cons(double *x, double *par){
        return par[0];
        //return 0.106947;
    }
    void constFit(){
        TF1 * fitfuc = new TF1("fitfuc", cons, 0,3,1);
        TH1F * Dataa = new TH1F("asdada","", 3, 0, 3);
        Dataa->SetBinContent(1, 0.12643);
        Dataa->SetBinContent(2, 0.117893);
        Dataa->SetBinContent(3, 0.106947);
        Dataa->SetBinError(1, 0.019732);
        Dataa->SetBinError(2, 0.0280342);
        Dataa->SetBinError(3, 0.0231889);
        fitfuc->FixParameter(0,0.106947);
        Dataa->Fit(fitfuc,"R");
    }