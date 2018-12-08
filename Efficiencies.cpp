#include<iostream>
#include<TMath.h>
#include"ErrorProp.hpp"
using namespace std;

#define EFFI_UP_POS  0.376192// +- 0.0400368
//!just coincident because when we change the shoulder position it changes drastically, 0.42 say.

#define EFFI_DOWN_POS  0.493563// +- 0.0154847
#define Err_Asymm 0.0229557


void Efficiencies()
{
    Double_t Effi_Pos = (EFFI_UP_POS + EFFI_DOWN_POS) / 2;
    Double_t Effi_Neg = 1-Effi_Pos;
    Double_t Asym_Acc = Effi_Neg-Effi_Pos;

    
    cout<<"\n\nEfficiency N+ =\t"<<Effi_Pos<<endl;
    cout<<"Efficiency N- =\t"<<Effi_Neg<<endl;
    cout<<"Effi Glb Asym =\t"<<Asym_Acc<<" +- "<<Err_Asymm<<endl;
}