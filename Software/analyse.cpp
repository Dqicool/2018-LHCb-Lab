#include "Analysis.hpp"

// This is the analysis class, which realises the generic Analysis
// from Analysis.hpp
//
// Look in Analysis.hpp for the event variables available.
class MyAnalysis : public Analysis {
public:
    // Define your histograms here
    TH1F        *h_PX;
    TH1F        *h_PY;
    TH1F        *h_PZ;
    TH2F        *h_TXTY;

    TH1F        *h_H1Pi;
    TH1F	    *h_H1Ka;
    TH1F	    *h_H2Pi;
    TH1F	    *h_H2Ka;
    TH1F	    *h_H3Pi;
    TH1F	    *h_H3Ka;

    void     BookHistos();

    Bool_t   Cut();
    void     Execute();
};

void MyAnalysis::BookHistos()
{
    // This function is only called once at the start of the program.
    // Book your histograms here. The format is object_name,
    // histogram_name, number_of_bins, minimum, maximum For a 2D
    // histogram, use TH2F with first the number of bins and limits
    // for the x axis and then for the y axis
    //
    // push_back() adds the histograms to a vector v_Histos.  This
    // will take care of writing out histograms in
    // Analysis::SaveHistos
    v_Histos.push_back( h_PX   = new TH1F("h_PX",  "", 100, -1e4, 1e4) );
    v_Histos.push_back( h_PY   = new TH1F("h_PY",  "", 100, -1e4, 1e4) );
    v_Histos.push_back( h_PZ   = new TH1F("h_PZ",  "", 100, -1e5, 1e5) );
    v_Histos.push_back( h_TXTY = new TH2F("h_TXTY","", 100, -1,1, 100,-1, 1) );

    v_Histos.push_back(h_H1Pi = new TH1F("h_H1Pi","", 100, 0, 1) );
    v_Histos.push_back(h_H1Ka = new TH1F("h_H1Ka","", 100, 0, 1) );
    v_Histos.push_back(h_H2Pi = new TH1F("h_H2Pi","", 100, 0, 1) );
    v_Histos.push_back(h_H2Ka = new TH1F("h_H2Ka","", 100, 0, 1) );
    v_Histos.push_back(h_H3Pi = new TH1F("h_H3Pi","", 100, 0, 1) );
    v_Histos.push_back(h_H3Ka = new TH1F("h_H3Ka","", 100, 0, 1) );
}

Bool_t MyAnalysis::Cut()
{
    // This function is called for every event from the Execute
    // function to define whether or not to accept this event.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    // This example checks if the PZ component of particle 3 is greater than 0. 


    if ( H3_PZ < 0 )
	return false;
    else
	return true;
}

void MyAnalysis::Execute()
{
    // This method gets called on every event.
    // In this example the momentum components are filled into histograms.

    // Call the Cut function to decide whether to plot this event or not
    // it returns if the cut function returns false
    if ( !Cut() )
	return;

    // Fill your histograms below.
    // fill the momentum of all three particles 
    h_PX->Fill( H1_PX );
    h_PX->Fill( H2_PX );
    h_PX->Fill( H3_PX );
    // the PY of all three particles
    h_PY->Fill( H1_PY );
    h_PY->Fill( H2_PY );
    h_PY->Fill( H3_PY );
    // the PZ of all three particles
    h_PZ->Fill( H1_PZ );
    h_PZ->Fill( H2_PZ );
    h_PZ->Fill( H3_PZ );
    // 2D histogram of PX/PZ vs PY/PZ
    h_TXTY->Fill( H1_PX / H1_PZ, H1_PY / H1_PZ );
    
    // Do same for probabilityies of being a pion/kaon
    h_H1Pi->Fill( H1_ProbPi );
    h_H1Ka->Fill( H1_ProbK );
    h_H2Pi->Fill( H2_ProbPi );
    h_H2Ka->Fill( H2_ProbK );
    h_H3Pi->Fill( H3_ProbPi );
    h_H3Ka->Fill( H3_ProbK );
}


// The main function just calls the generic AnalysisMain function
// with the MyAnalysis class
//
// Normally you don't need to change this
int main(int argc, char* argv[])
{
    MyAnalysis* ana = new MyAnalysis();
    int res = ana->AnalysisMain(argc, argv);
    return res;
}
