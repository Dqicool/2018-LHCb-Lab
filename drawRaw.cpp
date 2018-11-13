#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

void drawRaw(){
		TFile *f1 = new TFile("Output/DataAll.root");
		TFile *f2 = new TFile("Data/B2HHH_MagnetUp.root");
        TH1F *a = (TH1F*)f1->Get("a");
        TH1F *b = (TH1F*)f1->Get("b");
        TH1F *c = (TH1F*)f1->Get("c");
        TH1F *d = (TH1F*)f1->Get("d");
        TH1F *e = (TH1F*)f1->Get("e");
        //a->SetAxisRange(0,200);
        //b->SetAxisRange(0,12);
        c->SetAxisRange(0,50);
        d->SetAxisRange(0,50);
        e->SetAxisRange(0,50);
        a->SetAxisRange(0,50);
        b->SetAxisRange(0,20);

        TCanvas *c1= new TCanvas("c1","",600,600);
		d->SetLineColor(kBlue);
	    d->Draw();
		c->SetLineColor(kRed);
		c->Draw("same");
        e->SetLineColor(kGreen);
        e->Draw("same");

        TCanvas *c2= new TCanvas("c2","",600,600);
		a->Draw();
        
        TCanvas *c3= new TCanvas("c3","",600,600);
		b->Draw();
        
}