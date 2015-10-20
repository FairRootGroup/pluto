////////////////////////////////////////////////////////
//  PPropagator Class implementation file
//
//  Can be used for Pion-Nukleon-reactions (see ref. 23)
//
//                    Author:  H. Schuldes / IF
//                    Written: 21.04.09
//
////////////////////////////////////////////////////////

#include "PPropagator.h"


PPropagator::PPropagator(const Char_t *id, const Char_t *de, Int_t key) : 
    PChannelModel(id, de,key) {
    //Constructor
    pid = is_pid;
};


PDistribution* PPropagator::Clone(const char*) const
{
    //clone the object
    return new PPropagator((const PPropagator &)* this);
};

TComplex  PPropagator::GetAmplitude(Double_t *mass, Int_t *didx) {
    
    Double_t x = mass[0];
    
    Double_t M = makeStaticData()->GetParticleMass(pid); 
    Double_t Gamma = makeDynamicData()->GetParticleTotalWidth(mass[0],pid); 

//     if (pid == 52)
// 	Gamma = makeStaticData()->GetParticleTotalWidth(pid); 
//     else  {

// 	Gamma =  makeStaticData()->GetParticleTotalWidth(pid)
// 	    * (M/x) * sqrt(pow((x*x - 4*0.134*0.134) / (M*M-4*0.134*0.134)  ,3));
//     }
	

    Double_t Re = x*x - M*M;
    Double_t Im = Gamma * M;
    TComplex Denom(Re, Im); 
    TComplex Nom(1.,0);
    TComplex S = Nom/Denom; 
    
    return S;    
};

ClassImp(PPropagator)
