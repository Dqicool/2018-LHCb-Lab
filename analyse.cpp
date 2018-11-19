
#include "Analysis.hpp"
#define C    1	         //Natural units for c
#define M0K  493.667     //Mev
#define M0Pi 139.570     //MeV
#define M0D  1864.84     //MeV
#define M0B  5279.29     //Mev

#define HEI 4 //0.25 0.5 1 2 4
#define WID 4 //4


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
    TH1F *a,*b,*c,*d,*e;

    TH1F        *h_H1Pi;
    TH1F	    *h_H1Ka;
    TH1F	    *h_H2Pi;
    TH1F	    *h_H2Ka;
    TH1F	    *h_H3Pi;
    TH1F	    *h_H3Ka;
    
    TH1F	    *h_B_M0_Pos;
    TH1F        *h_B_M0_Neg;
    
    TH1F	    *h_H_M;
    TH1F	    *h_L_M;
    
    TH2F	    *h_Dalitz_Pos_Com;
    TH2F        *h_Dalitz_Pos_Bac;
    TH2F	    *h_Dalitz_Neg_Com;
    TH2F        *h_Dalitz_Neg_Bac;


    void     BookHistos();

    //Bool_t   Cut();
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
    v_Histos.push_back( h_PX   = new TH1F("h_PX",  "", 1000, -1e4, 1e4) );
    v_Histos.push_back( h_PY   = new TH1F("h_PY",  "", 1000, -1e4, 1e4) );
    v_Histos.push_back( h_PZ   = new TH1F("h_PZ",  "", 1000, -1e5, 1e5) );
    v_Histos.push_back( h_TXTY = new TH2F("h_TXTY","", 1000, -1,1, 1000,-1, 1) );
    v_Histos.push_back( a = new TH1F("B_FlightDistance", "",1000,0,100));
    v_Histos.push_back( b = new TH1F("B_VertexChi2", "",1000,0,100));
    v_Histos.push_back( c = new TH1F("H1_IPChi2", "",1000,0,100));
    v_Histos.push_back( d = new TH1F("H2_IPChi2", "",1000,0,100));
    v_Histos.push_back( e = new TH1F("H3_IPChi2", "",1000,0,100));

    //Pion Kaon Possibility
    v_Histos.push_back(h_H1Pi = new TH1F("h_H1Pi","", 1000, 0, 1) );
    v_Histos.push_back(h_H1Ka = new TH1F("h_H1Ka","", 1000, 0, 1) );
    v_Histos.push_back(h_H2Pi = new TH1F("h_H2Pi","", 1000, 0, 1) );
    v_Histos.push_back(h_H2Ka = new TH1F("h_H2Ka","", 1000, 0, 1) );
    v_Histos.push_back(h_H3Pi = new TH1F("h_H3Pi","", 1000, 0, 1) );
    v_Histos.push_back(h_H3Ka = new TH1F("h_H3Ka","", 1000, 0, 1) );
    //Invariant Mass of B
    v_Histos.push_back(h_B_M0_Pos =  new TH1F("h_B_M0_Pos",  "", 250, 5e3, 6e3) );
    v_Histos.push_back(h_B_M0_Neg =  new TH1F("h_B_M0_Neg",  "", 250, 5e3, 6e3) );
    //Higher and Lower Mass Pairs
    v_Histos.push_back(h_H_M =  new TH1F("h_H_M",  "", 500, 0, 5.8e3) );
    v_Histos.push_back(h_L_M =  new TH1F("h_L_M",  "", 500, 0, 5.8e3) );

    //Dalitz mass data
    v_Histos.push_back(h_Dalitz_Pos_Com = new TH2F("h_Dalitz_Pos_Com", "", 24/WID, -4, 20, 32/HEI, -2, 30));
    v_Histos.push_back(h_Dalitz_Pos_Bac = new TH2F("h_Dalitz_Pos_Bac", "", 24/WID, -4, 20, 32/HEI, -2, 30));
    v_Histos.push_back(h_Dalitz_Neg_Com = new TH2F("h_Dalitz_Neg_Com", "", 24/WID, -4, 20, 32/HEI, -2, 30));
    v_Histos.push_back(h_Dalitz_Neg_Bac = new TH2F("h_Dalitz_Neg_Bac", "", 24/WID, -4, 20, 32/HEI, -2, 30));

}

/*Bool_t MyAnalysis::Cut() // cut() has been abandoned
{
    // This function is called for every event from the Execute
    // function to define whether or not to accept this event.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    // This example checks if the PZ component of particle 3 is greater than 0. 


    return true;
}*/

void MyAnalysis::Execute(){
    // This method gets called on every event.
    // In this example the momentum components are filled into histograms.


    //equal functionality to cut() 
        if ( (H1_ProbPi < 0.843) || (H2_ProbPi < 0.843) || (H3_ProbPi < 0.843))
            return;
        if ( H1_isMuon || H2_isMuon || H3_isMuon)
            return;
        //if (H1_ProbPi<H1_ProbK || H2_ProbPi<H2_ProbK||H3_ProbPi<H3_ProbK)
        //    return;
    // fill abcde
        a->Fill(B_FlightDistance);
        b->Fill(B_VertexChi2);
        c->Fill(H1_IPChi2);
        d->Fill(H2_IPChi2);
        e->Fill(H3_IPChi2);
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
    
    // Do same for probabilities of being a pion/kaon
        h_H1Pi->Fill( H1_ProbPi );
        h_H1Ka->Fill( H1_ProbK );
        h_H2Pi->Fill( H2_ProbPi );
        h_H2Ka->Fill( H2_ProbK );
        h_H3Pi->Fill( H3_ProbPi );
        h_H3Ka->Fill( H3_ProbK );
    
    //Invariant Mass of B
        //Calculate total momentum for each decay particle
        float P1 = TMath::Sqrt( pow(H1_PX,2) + pow(H1_PY,2) + pow(H1_PZ,2));
        float P2 = TMath::Sqrt( pow(H2_PX,2) + pow(H2_PY,2) + pow(H2_PZ,2));
        float P3 = TMath::Sqrt( pow(H3_PX,2) + pow(H3_PY,2) + pow(H3_PZ,2));
        //Calculate total energy for each decay particle using known invariant mass
        float E1 = TMath::Sqrt(pow(P1*C,2) + pow(M0Pi*C*C,2));
        float E2 = TMath::Sqrt(pow(P2*C,2) + pow(M0Pi*C*C,2));
        float E3 = TMath::Sqrt(pow(P3*C,2) + pow(M0Pi*C*C,2));
        //Calculate the x, y, z components of the B-meson momentum 
        float P_BX = H1_PX+H2_PX+H3_PX;
        float P_BY = H1_PY+H2_PY+H3_PY;
        float P_BZ = H1_PZ+H2_PZ+H3_PZ;
        //Calculate the B-meson absolute momentum
        float abs_P_B = TMath::Sqrt( pow(P_BX,2) + pow(P_BY,2) + pow(P_BZ,2));
        //Calculate the invariant mass of the B-meson
        float M0_B = TMath::Sqrt(pow(E1+E2+E3,2)-pow(abs_P_B*C,2))/pow(C,2);
        //Invariant Mass of B 

    //High and Low Mass Pair
        float M1, M2, MH, ML;
        //select B+ Decay
        if (H1_Charge+H2_Charge+H3_Charge == 1)
            {
                //Choosing Mass Pair
                if (H1_Charge == -1){
                    M1 = TMath::Sqrt(pow(E1+E2,2)-(pow(H1_PX + H2_PX,2)+pow(H1_PY + H2_PY,2)+pow(H1_PZ + H2_PZ,2)));
                    M2 = TMath::Sqrt(pow(E1+E3,2)-(pow(H1_PX + H3_PX,2)+pow(H1_PY + H3_PY,2)+pow(H1_PZ + H3_PZ,2)));
                }
                else if (H2_Charge == -1){
                    M1 = TMath::Sqrt(pow(E2+E1,2)-(pow(H2_PX + H1_PX,2)+pow(H2_PY + H1_PY,2)+pow(H2_PZ + H1_PZ,2)));
                    M2 = TMath::Sqrt(pow(E2+E3,2)-(pow(H2_PX + H3_PX,2)+pow(H2_PY + H3_PY,2)+pow(H2_PZ + H3_PZ,2)));
                }
                else if (H3_Charge == -1){
                    M1 = TMath::Sqrt(pow(E3+E1,2)-(pow(H3_PX + H1_PX,2)+pow(H3_PY + H1_PY,2)+pow(H3_PZ + H1_PZ,2)));
                    M2 = TMath::Sqrt(pow(E3+E2,2)-(pow(H3_PX + H2_PX,2)+pow(H3_PY + H2_PY,2)+pow(H3_PZ + H2_PZ,2)));
                }
                if (M2 > M1){ MH = M2; ML = M1; }
                else        { MH = M1; ML = M2; }
                if (pow(M0_B-M0B,2) < 1369)
	                h_Dalitz_Pos_Com->Fill(ML*ML/1e6,MH*MH/1e6);
                if(M0_B >5400 && M0_B<5770)
                    h_Dalitz_Pos_Bac->Fill(ML*ML/1e6*(M0B/M0_B)*(M0B/M0_B),MH*MH/1e6*(M0B/M0_B)*(M0B/M0_B));
            }
            
        //select B- Decay
        else if (H1_Charge+H2_Charge+H3_Charge == -1)
            {
                //Choosing Mass Pair
                if (H1_Charge == 1){
                    M1 = TMath::Sqrt(pow(E1+E2,2)-(pow(H1_PX + H2_PX,2)+pow(H1_PY + H2_PY,2)+pow(H1_PZ + H2_PZ,2)));
                    M2 = TMath::Sqrt(pow(E1+E3,2)-(pow(H1_PX + H3_PX,2)+pow(H1_PY + H3_PY,2)+pow(H1_PZ + H3_PZ,2)));
                }
                else if (H2_Charge == 1){
                    M1 = TMath::Sqrt(pow(E2+E1,2)-(pow(H2_PX + H1_PX,2)+pow(H2_PY + H1_PY,2)+pow(H2_PZ + H1_PZ,2)));
                    M2 = TMath::Sqrt(pow(E2+E3,2)-(pow(H2_PX + H3_PX,2)+pow(H2_PY + H3_PY,2)+pow(H2_PZ + H3_PZ,2)));
                }
                else if (H3_Charge == 1){
                    M1 = TMath::Sqrt(pow(E3+E1,2)-(pow(H3_PX + H1_PX,2)+pow(H3_PY + H1_PY,2)+pow(H3_PZ + H1_PZ,2)));
                    M2 = TMath::Sqrt(pow(E3+E2,2)-(pow(H3_PX + H2_PX,2)+pow(H3_PY + H2_PY,2)+pow(H3_PZ + H2_PZ,2)));
                }
                if (M2 > M1){ MH = M2; ML = M1; }
                else        { MH = M1; ML = M2; }
                if (pow(M0_B-M0B,2) < 1369)
	                h_Dalitz_Neg_Com->Fill(ML*ML/1e6,MH*MH/1e6);
                if(M0_B >5400 && M0_B<5770)
                    h_Dalitz_Neg_Bac->Fill(ML*ML/1e6*(M0B/M0_B)*(M0B/M0_B),MH*MH/1e6*(M0B/M0_B)*(M0B/M0_B));
            }

        //Discard the D meson Decay by Simply cut M1 or M2 in D meson Mass
        if (pow(MH - M0D ,2) > 2500 && pow(ML-M0D, 2) > 2500){
            if(H1_Charge+H2_Charge+H3_Charge == 1){
                h_B_M0_Pos->Fill( M0_B );}
            else{
                h_B_M0_Neg->Fill( M0_B );}

            h_H_M->Fill( MH );
            h_L_M->Fill( ML );
        }

}


// The main function just calls the generic AnalysisMain function
// with the MyAnalysis class
//
// Normally you don't need to change this
int main(int argc, char* argv[]) {
    MyAnalysis* ana = new MyAnalysis();
    int res = ana->AnalysisMain(argc, argv);
    return res;
}
