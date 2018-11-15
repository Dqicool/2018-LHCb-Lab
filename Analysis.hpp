#ifndef _ANALYSIS_HPP
#define _ANALYSIS_HPP

#include <iostream>
#include <cstdlib>
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TF1.h"


class Analysis {
    public:
    // Declaration of leaf types
        Double_t         B_FlightDistance;
        Double_t         B_VertexChi2;
        Double_t         H1_PX;
        Double_t         H1_PY;
        Double_t         H1_PZ;
        Double_t         H1_ProbK;
        Double_t         H1_ProbPi;
        Int_t            H1_Charge;
        Double_t         H1_IPChi2;
        Int_t            H1_isMuon;
        Double_t         H2_PX;
        Double_t         H2_PY;
        Double_t         H2_PZ;
        Double_t         H2_ProbK;
        Double_t         H2_ProbPi;
        Int_t            H2_Charge;
        Double_t         H2_IPChi2;
        Int_t            H2_isMuon;
        Double_t         H3_PX;
        Double_t         H3_PY;
        Double_t         H3_PZ;
        Double_t         H3_ProbK;
        Double_t         H3_ProbPi;
        Int_t            H3_Charge;
        Double_t         H3_IPChi2;
        Int_t            H3_isMuon;

    // tree and vector of histos
        TChain         *myChain;
        std::vector< TNamed* > v_Histos;

        void     Init(TChain *chain, std::string choice);
        Int_t    GetEntry(Long64_t entry);
        void     Loop(int nevts = -1);
        void     SaveHistos(const char*);

        int      AnalysisMain(int argc, char* argv[]);    

        // These are pure virtual -- should be defined in real analysis child class
        virtual void     BookHistos() = 0;
        virtual void     Execute() = 0;
};

Int_t Analysis::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!myChain) return 0;
    return myChain->GetEntry(entry);
}

void Analysis::SaveHistos(const char* fname)
{
    TFile *f = new TFile(fname, "RECREATE");
    f->cd();
    std::vector< TNamed* >::iterator it = v_Histos.begin();
    for( ; it != v_Histos.end(); it++ ) {
	(*it)->Write();
    }  
    f->Close();
}

void Analysis::Init(TChain *chain, std::string choice)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!chain) return;
    myChain = chain;
    myChain->SetMakeClass(1);

    myChain->SetBranchAddress("H1_PX", &H1_PX);
    myChain->SetBranchAddress("H1_PY", &H1_PY);
    myChain->SetBranchAddress("H1_PZ", &H1_PZ);
    myChain->SetBranchAddress("H2_PX", &H2_PX);
    myChain->SetBranchAddress("H2_PY", &H2_PY);
    myChain->SetBranchAddress("H2_PZ", &H2_PZ);
    myChain->SetBranchAddress("H3_PX", &H3_PX);
    myChain->SetBranchAddress("H3_PY", &H3_PY);
    myChain->SetBranchAddress("H3_PZ", &H3_PZ);
    
    if ( "PhaseSpace" == choice ) {
	// some variables don't exist, so set them to default values
	    B_FlightDistance = 0.;
        B_VertexChi2 = 0.;
        H1_ProbK = 0.;
        H1_ProbPi = 1;
        H1_Charge = 1;
        H1_IPChi2 = 0.;
        H1_isMuon = 0;
        H2_ProbK = 0.;
        H2_ProbPi = 1.;
        H2_Charge = 1;
        H2_IPChi2 = 0.;
        H2_isMuon = 0;
        H3_ProbK = 0.;
        H3_ProbPi = 1.;
        H3_Charge = -1;
        H3_IPChi2 = 0.;
        H3_isMuon = 0;
    }
    else {
	// normal tree, declare remaining branches
        myChain->SetBranchAddress("B_FlightDistance", &B_FlightDistance);
        myChain->SetBranchAddress("B_VertexChi2", &B_VertexChi2);
        myChain->SetBranchAddress("H1_ProbK", &H1_ProbK);
        myChain->SetBranchAddress("H1_ProbPi", &H1_ProbPi);
        myChain->SetBranchAddress("H1_Charge", &H1_Charge);
        myChain->SetBranchAddress("H1_IPChi2", &H1_IPChi2);
        myChain->SetBranchAddress("H1_isMuon", &H1_isMuon);
        myChain->SetBranchAddress("H2_ProbK", &H2_ProbK);
        myChain->SetBranchAddress("H2_ProbPi", &H2_ProbPi);
        myChain->SetBranchAddress("H2_Charge", &H2_Charge);
        myChain->SetBranchAddress("H2_IPChi2", &H2_IPChi2);
        myChain->SetBranchAddress("H2_isMuon", &H2_isMuon);
        myChain->SetBranchAddress("H3_ProbK", &H3_ProbK);
        myChain->SetBranchAddress("H3_ProbPi", &H3_ProbPi);
        myChain->SetBranchAddress("H3_Charge", &H3_Charge);
        myChain->SetBranchAddress("H3_IPChi2", &H3_IPChi2);
        myChain->SetBranchAddress("H3_isMuon", &H3_isMuon);
    }
    BookHistos();
}

void Analysis::Loop(int nevts)
{
    //   In a ROOT session, you can do:
    //      Root > .L Analysis.C
    //      Root > Analysis t
    //      Root > t.GetEntry(12); // Fill t data members with entry number 12
    //      Root > t.Loop();       // Loop on all entries
    //

    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry
    //
    //       To read only selected branches, Insert statements like:
    // METHOD1:
    //    myTree->SetBranchStatus("*",0);  // disable all branches
    //    myTree->SetBranchStatus("branchname",1);  // activate branchname
    // METHOD2: replace line
    //    myTree->GetEntry(jentry);       //read all branches
    //by  b_branchname->GetEntry(ientry); //read only this branch
    if (myChain == 0) return;

    Long64_t nentries = myChain->GetEntries();
    // limit number of events to those specified on the command prompt. Default is to read all.
    if ( ( nevts > 0 ) && ( nevts < nentries ) ) nentries = nevts;
    std::cout << "Performing the analysis on " << nentries << " events" << std::endl;

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t ientry=0; ientry<nentries;ientry++) {
	if (( 0 == ientry % 100000 ) || (ientry==nentries-1)) std::cout << "At entry " << ientry << std::endl;
	nb = myChain->GetEntry(ientry);   nbytes += nb;
	Execute();
    }
}

int Analysis::AnalysisMain(int argc, char* argv[])
{
    // setting analysis type and configuring input data
    std::string oname = "";
    std::string choice = "";
    TChain* chain;
    if ( argc > 1 ) {
        choice = argv[1];
        if ( "PhaseSpace" == choice ) { 
            chain = new TChain("PhaseSpaceTree");
            chain->Add( "Data/PhaseSpaceSimulation.root" );
            oname = "Output/PhaseSpace.root";
            std::cout << "Performing phase space analysis" << std::endl;
        }    
        else if ( "DataMagnetDown" == choice ) { 
            chain = new TChain("DecayTree");
            chain->Add( "Data/B2HHH_MagnetDown.root" );
            oname = "Output/DataMagnetDown.root";
            std::cout << "Performing data analysis with magnet polarity down" << std::endl;
        }    
        else if ( "DataMagnetUp" == choice ) { 
            chain = new TChain("DecayTree");
            chain->Add( "Data/B2HHH_MagnetUp.root" );
            oname = "Output/DataMagnetUp.root";
            std::cout << "Performing data analysis with magnet polarity up" << std::endl;
        }    
        else if ( "DataAll" == choice ) { 
            chain = new TChain("DecayTree");
            chain->Add( "Data/B2HHH_MagnetDown.root" );
            chain->Add( "Data/B2HHH_MagnetUp.root" );
            oname = "Output/DataAll.root";
            std::cout << "Performing data analysis with both magnet polarities" << std::endl;
        }    
        else {
            std::cout << "Unknown analysis type" << std::endl;
            std::cout << "Options are PhaseSpace, DataAll, DataMagnetDown, DataMagnetUp" << std::endl;
            return 0;
        }
    }
    else {
        std::cout << "Please specify which samples you wish to analyse" << std::endl;
        std::cout << "Options are PhaseSpace, DataAll, DataMagnetDown, DataMagnetUp" << std::endl;
        std::cout << "You many also specify the number of events to analyse. If this argument is not given the full sample is analysed." << std::endl;
        return 0;
    }

    // checking whether number of events is to be limited
    int nevts = -1;
    if ( argc > 2 ) { 
        nevts = atoi( argv[2] );
        std::cout << "Setting number of output events to " << nevts << std::endl;
    }
    Init(chain, choice);
    Loop(nevts);
    SaveHistos(oname.c_str());

    return 0;
}

#endif