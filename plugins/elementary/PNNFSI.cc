/////////////////////////////////////////////////////////////////////
//
// NN FSI from the Jost function
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PNNFSI.h"
#include "PDynamicData.h"


ClassImp(PNNFSI)
    
    PNNFSI::PNNFSI()  {

    Setup();

} ;

PNNFSI::PNNFSI(const Char_t *id,const  Char_t *de, Int_t key) :
    PChannelModel(id, de, key) {

    Setup();

} ;

void PNNFSI::Setup(void) {
    n1 = n2 = NULL;
    daughter_pos=0;
    for (int i=0; i< NNFSI_MAX_DAUGHTERS;i++) {
	daughters[i]= NULL;
    }

    triplet =0;
    a0=-23.768;
    r0=2.75;
    is_np =0;

}

PDistribution* PNNFSI::Clone(const char*delme) const {
    return new PNNFSI((const PNNFSI &)* this);
};


Bool_t PNNFSI::Init(void) {

    //looking for 2 nucleons. This is mandatory
    n1 = GetParticle("p");
    n2 = GetParticle("p");
    if (!n2) {
	n2 = GetParticle("n");
	is_np =1;
    }

    if (!n1 || !n2) {
	Warning("Init","No 2 nucleons");
	return kFALSE;
    }

    //Now get all daugthers
    for (int i=0; i<NNFSI_MAX_DAUGHTERS;i++) {
	daughters[i]= GetParticle("daughter");
	if (daughters[i]) {	    
	    daughter_pos++;
	}
    }

    return kTRUE;    
};

Double_t PNNFSI::GetWeight(Double_t *mass, Int_t *) {
    TComplex jost(1.,0.);
    Double_t k= mass[0]/0.197 ;
    Double_t alpha = 1./r0 * (
	sqrt (1.-2.*r0/a0) +1.
	);
    Double_t beta = 1./r0 * (
	sqrt (1.-2.*r0/a0) -1.
	);
    if (!triplet) {
	TComplex num(k,beta);
	TComplex denom(k,alpha);
	jost = num/denom;
    } else {
	cout << "TODO" << endl;
    }

//    if (!is_np)
    return (1./jost.Rho2());

}

Double_t PNNFSI::GetWeight(void) {
//     n1->Print();
//     n2->Print();

    TLorentzVector a = ((TLorentzVector) *n1) + ((TLorentzVector) *n2);
    TLorentzVector b = ((TLorentzVector) *n1);

    b.Boost(-a.BoostVector());

    Double_t k=b.P();
    weight = GetWeight(&k);
    // cout << k << ":" << weight << endl;
    return weight;


}



ClassImp(PNNFSI)

