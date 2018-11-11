#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"

double_t Background(double_t *x, double_t *par)
{
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

double_t LorentzPeak(double_t *x, double_t *par)
{
    return (0.5*par[0]*par[1]/TMath::Pi()) / TMath::Max(1.e-10, (x[0]-par[2])*(x[0]-par[2])+ 0.25*par[1]*par[1]);
}

double_t SumOfFunc(double_t *x, double_t *par)
{ 
    double_t val = Background(x,&par[0]) + LorentzPeak(x, &par[3]);
    return val;
}

void ComFun()
{
    //construct histo
        const int nBins = 60;
        Stat_t data[nBins] = {   6, 1,10,12, 6,13,23,22,15,21, 
                                23,26,36,25,27,35,40,44,66,81, 
                                75,57,48,45,46,41,35,36,53,32, 
                                40,37,38,31,36,44,42,37,32,32, 
                                43,44,35,33,33,39,29,41,32,44, 
                                26,39,29,35,32,21,21,15,25,15 };
        TH1F *histo = new TH1F("example_9_1", "Lorentzian Peak on Quadratic Background",60,0,3);
        for(int i=0; i < nBins; i++) {
            // we use these methods to explicitly set the content 
            // and error instead of using the fill method. 
            histo->SetBinContent(i+1,data[i]); 
            histo->SetBinError(i+1,TMath::Sqrt(data[i]));
        }
    //Create fit function
        TF1 *fitfunc = new TF1("fitfunc",SumOfFunc, 0,3,6);
    //fitting
        fitfunc->SetParameters(1,1,1,1,1,1);
        histo->Fit("fitfunc");
        double_t par[6];
        fitfunc->GetParameters(par);
        TF1 *bacfun = new TF1("bacfun",Background,0,3,3);
        TF1 *lorfun = new TF1("lorfun",LorentzPeak,0,3,3);
        bacfun->SetParameters(&par[0]);
        lorfun->SetParameters(&par[3]);
        lorfun->Draw("same");
        bacfun->Draw("same");
}



