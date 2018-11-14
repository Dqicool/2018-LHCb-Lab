    #include <TFile.h>
    #include <TH1.h>
    #include <TH2.h>
    #include <TCanvas.h>
    #include <TMath.h>
    #include <TF1.h>
    #include <iostream>
    using namespace std;

Double_t func(Double_t *x, Double_t *par){
    if ( x[0]<= par[0])
        return par[5]*TMath::Exp(-(pow(x[0]-par[0],2)/(2*(pow(par[1],2)+par[2]*(pow(x[0]-par[0],2))))));
    else
        return par[5]*TMath::Exp(-(pow(x[0]-par[0],2)/(2*(pow(par[3],2)+par[4]*(pow(x[0]-par[0],2))))));
}
void cruijffVar()
{
    TF1 *fun1 = new TF1("fun1", func, -10,10,6);
    fun1->SetParNames( "x0", "SigmaL", "AlphaL","SigmaR", "AlphaR", "I");
    fun1->SetParameters(0,    2,        0.1,       2,        0.2,    10);
    fun1->Draw();
}