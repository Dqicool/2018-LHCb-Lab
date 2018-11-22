#include<iostream>
#include<TMath.h>
using namespace std;

#define EFFI_UP_POS  0.4444943
#define ERR_EFFI_UP_POS 0.0647
#define EFFI_DOWN_POS  0.444798
#define ERR_EFFI_DOWN_POS 0.0710
#define EFFI_UP_NEG  1 - EFFI_UP_POS
#define ERR_EFFI_UP_NEG 0.0812351
#define EFFI_DOWN_NEG  1 - EFFI_DOWN_POS
#define ERR_EFFI_DOWN_NEG 0.0835704

// Error Propagation functions
    Double_t ErrAPlusB(Double_t ErrA, Double_t ErrB)
    {
        Double_t Err2 = pow(ErrA,2)+pow(ErrB,2); 
        return TMath::Sqrt(Err2);
    }

void Efficiencies()
{
    Double_t Effi_Pos = (EFFI_UP_POS + EFFI_DOWN_POS) / 2;
    Double_t Effi_Neg = (EFFI_UP_NEG + EFFI_DOWN_NEG) / 2;
    Double_t Asym_Acc = Effi_Neg-Effi_Pos;
    Double_t Err_Effi_Pos = ErrAPlusB(ERR_EFFI_UP_POS, ERR_EFFI_DOWN_POS)/2;
    Double_t Err_Effi_Neg = ErrAPlusB(ERR_EFFI_UP_NEG, ERR_EFFI_DOWN_NEG)/2;
    Double_t Err_Asym_Acc = ErrAPlusB(Err_Effi_Pos,Err_Effi_Neg);
    
    cout<<"Efficiency N+ =\t"<<Effi_Pos<<" +- "<<Err_Effi_Pos<<endl;
    cout<<"Efficiency N- =\t"<<Effi_Neg<<" +- "<<Err_Effi_Neg<<endl;
    cout<<"Effi Glb Asym =\t"<<Asym_Acc<<" +- "<<Err_Asym_Acc<<endl;
}