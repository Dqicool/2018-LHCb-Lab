#include <iostream>
#include <stdlib.h>
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TF1.h"

#ifndef _ERROR_PROP_
#define _ERROR_PROP_
//6 Error Propagation functions
    Double_t ErrAPlusB(Double_t ErrA, Double_t ErrB)
    {
        Double_t Err2 = pow(ErrA,2)+pow(ErrB,2)+2*ErrA*ErrB; 
        return TMath::Sqrt(Err2);
    }
    Double_t ErrAMinusB(Double_t ErrA, Double_t ErrB)
    {
        Double_t Err2 = pow(ErrA,2) + pow(ErrB,2) - 2*ErrA*ErrB; 
        return TMath::Sqrt(Err2);
    }
    Double_t ErrAMultB(Double_t C, Double_t A, Double_t B, Double_t ErrA, Double_t ErrB)
    {
        Double_t f2 = C*C;
        Double_t AErrA2 = pow(ErrA/A,2);
        Double_t BErrB2 = pow(ErrB/B,2);
        Double_t CoriAB = ErrA*ErrB/(B*A);
        Double_t Err2 = f2*(AErrA2+BErrB2+2*CoriAB);
        return TMath::Sqrt(Err2);
    }
    Double_t ErrADiviB(Double_t C, Double_t A, Double_t B, Double_t ErrA, Double_t ErrB)
    {
        Double_t f2 = C*C;
        Double_t AErrA2 = pow(ErrA/A,2);
        Double_t BErrB2 = pow(ErrB/B,2);
        Double_t CoriAB = ErrA/A*ErrB/B;
        Double_t Err2 = f2*(AErrA2+BErrB2-2*CoriAB);
        return TMath::Sqrt(Err2);
    }
#endif