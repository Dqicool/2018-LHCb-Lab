#include "Analysis.hpp"

void drawOutput() {
    //////////////////////
    // Example Root Macro for third year B->hhh Lab
    ////////////////////// 

    // Open the root file that was produced by running the example program
    //TFile *f = new TFile("outputPhaseSpace.root");
    TFile *f1 = new TFile("Output/DataMagnetDown.root");
    TFile *f2 = new TFile("Output/PhaseSpace.root");
    TFile *f3 = new TFile("Output/DataAll.root");
    TFile *f4 = new TFile("Output/DataMagnetDown.root");

    // Get pointers to the example histograms that were made 
    TH1F *hx = (TH1F*)f1->Get("h_PX");
    TH1F *hy = (TH1F*)f1->Get("h_PY");
    TH2F *hxy = (TH2F*)f1->Get("h_TXTY");
    
    TH1F *H1Pi = (TH1F*)f1->Get("h_H1Pi");
    TH1F *H1Ka = (TH1F*)f1->Get("h_H1Ka");
    TH1F *H2Pi = (TH1F*)f1->Get("h_H2Pi");
    TH1F *H2Ka = (TH1F*)f1->Get("h_H2Ka");
    TH1F *H3Pi = (TH1F*)f1->Get("h_H3Pi");
    TH1F *H3Ka = (TH1F*)f1->Get("h_H3Ka");
    
    TH1F *hBM0 = (TH1F*)f3->Get("h_B_M0");

    TH1F *hHM  = (TH1F*)f3->Get("h_H_M");
    TH1F *hLM  = (TH1F*)f3->Get("h_L_M");
    
    // Create a canvas onto which the histograms are plotted and which can be saved
    TCanvas *c1 = new TCanvas("c1","",600,400);
    // Draw the first histogram with a blue line, and an x-axis title
    hx->SetLineColor(kBlue);
    hx->GetXaxis()->SetTitle("Momentum [MeV/c^{2}]");
    hx->Draw();
    // Draw the second histogram on the same plot with a red line
    hy->SetLineColor(kRed);
    hy->Draw("same");
    // Save the canvas as pdf file. Other possible formats include root (for later editing) eps, png.
    c1->SaveAs("Plots/hx_hy.pdf");

    // Repeat the above for the 2D histogram
      TCanvas *c2 = new TCanvas("c2","",600,400);
      hxy->SetStats(0);                           // remove the statistics box
      hxy->GetXaxis()->SetTitle("Slope in x");    // add axis titles
      hxy->GetYaxis()->SetTitle("Slope in y");
      hxy->Draw("colz");                          // draw with a colour scale
      c2->SaveAs("Plots/hxy.pdf");
      
      
    // Repeat for particle 1 probabilities
    TCanvas *c3 = new TCanvas("c3","",600,400);
    H1Pi->SetLineColor(kBlue);
    H1Pi->GetXaxis()->SetTitle("Probability of Particle 1 State");
    H1Pi->Draw();
    H1Ka->SetLineColor(kRed);
    H1Ka->Draw("same");
    c3->SaveAs("Plots/Particle1Prob.pdf");
    
    // Repeat for particle 2  
    TCanvas *c4 = new TCanvas("c4","",600,400);
    H2Pi->SetLineColor(kBlue);
    H2Pi->GetXaxis()->SetTitle("Probability of Particle 2 State");
    H2Pi->Draw();
    H2Ka->SetLineColor(kRed);
    H2Ka->Draw("same");
    c4->SaveAs("Plots/Particle2Prob.pdf");
    
    // Repeat for particle 3
    TCanvas *c5= new TCanvas("c5","",600,400);
    H3Pi->SetLineColor(kBlue);
    H3Pi->GetXaxis()->SetTitle("Probability of Particle 3 State");
    H3Pi->Draw();
    H3Ka->SetLineColor(kRed);
    H3Ka->Draw("same");
    c5->SaveAs("Plots/Particle3Prob.pdf");
    
    //Save plot of invariant mass for B meson
    TCanvas *c6= new TCanvas("c6","",600,400);
    hBM0->SetLineColor(kBlue);
    hBM0->GetXaxis()->SetTitle("B-Meson Invariant Mass");
    hBM0->Draw();
    c6->SaveAs("Plots/B-MesonInvariantMass_final_pion.pdf");

    TCanvas *c7= new TCanvas("c7","",600,400);
    hHM->SetLineColor(kBlue);
    hHM->GetXaxis()->SetTitle("B-Meson Invariant Mass");
    hHM->Draw();
    hLM->SetLineColor(kRed);
    hLM->Draw("same");
    c7->SaveAs("Plots/Diff_Mass_Pion_Pairs.pdf");
}
