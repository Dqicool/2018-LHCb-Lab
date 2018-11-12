#include"Analysis.hpp"

Double_t fitf(Double_t *x,Double_t *par) 
{
    Double_t arg = 0;
    if (par[2]!=0) arg = (x[0] - par[1])/par[2];
    Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg); 
    return fitval;
}


void SingleFit()
{
    TFile *f = new TFile("Output/DataAll.root"); 
    TH1F *hHM  = (TH1F*)f->Get("h_H_M");
    // Create a TF1 object using the function defined above.
    // The last three parameters specify the number of parameters // for the function.
    TF1 *func = new TF1("fit",fitf,-3,3,3);
    // set the parameters to the mean and RMS of the histogram 
    func->SetParameters(500,hHM->GetMean(),hHM->GetRMS());
    // give the parameters meaningful names
    func->SetParNames ("Constant","Mean","Sigma");
    // call TH1::Fit with the name of the TF1 object
    hHM->Fit("fit","","",3800,5000);
}

void Multifit()
{
    const Int_t np = 49;
    Float_t x[np] = {1.913521, 1.953769, 2.347435, 2.883654, 3.493567,
                    4.047560, 4.337210, 4.364347, 4.563004, 5.054247,
                    5.194183, 5.380521, 5.303213, 5.384578, 5.563983,
                    5.728500, 5.685752, 5.080029, 4.251809, 3.372246,
                    2.207432, 1.227541, 0.8597788,0.8220503,0.8046592,
                    0.7684097,0.7469761,0.8019787,0.8362375,0.8744895,
                    0.9143721,0.9462768,0.9285364,0.8954604,0.8410891,
                    0.7853871,0.7100883,0.6938808,0.7363682,0.7032954,
                    0.6029015,0.5600163,0.7477068,1.188785, 1.938228,
                    2.602717, 3.472962, 4.465014, 5.177035};
    TF1 *g1 = new TF1("m1","gaus",85,95);
    TF1 *g2 = new TF1("m2","gaus",98,108);
    TF1 *g3 = new TF1("m3","gaus",110,121);
    // The total is the sum of the three, each has 3 parameters 
    TF1 *total = new TF1("mstotal","gaus(0)+gaus(3)+gaus(6)",85,125);

    TH1F *h = new TH1F("g1","Example of several fits in subranges", np,85,134);
    h->SetMaximum(7);
    for (int i=0; i<np; i++) {
        h->SetBinContent(i+1,x[i]);
    }
    // Define the parameter array for the total function
    Double_t par[9];

    h->Fit(g1, "R");
    h->Fit(g2, "R+");
    h->Fit(g3, "R+");

    // Get the parameters from the fit
    g1->GetParameters(&par[0]);
    g2->GetParameters(&par[3]);
    g3->GetParameters(&par[6]);
    // Use the parameters on the sum
    total->SetParameters(par);
    h->Fit(total,"R+");
}

void CombinFunc()
{

}
void LearnFit()
{
    SingleFit();
    Multifit();
}