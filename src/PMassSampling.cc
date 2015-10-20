/////////////////////////////////////////////////////////////////////
//
// General purpose resonance mass sampling 
// Resonance shapes can be defined by a TF1-object
// This class allows to "play" around with different resonance shapes
// in coupled-channel calculations
//
//                                  Author:  Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PMassSampling.h"


ClassImp(PMassSampling)

PMassSampling::PMassSampling()  {

} ;

PMassSampling::PMassSampling(Char_t *id, Char_t *de, Int_t key) :
    PHadronModel(id, de, key) {
    
    shape1=NULL;

};

PDistribution* PMassSampling::Clone(const char*delme) const {
    return new PMassSampling((const PMassSampling &)* this);
};

Bool_t PMassSampling::SampleMass(Double_t *mass, Int_t *didx) {
    if (!shape1) {
	Warning("SampleMass","(%s): no TF1 found",GetDescription());
	mass[0] = 0;
	return kFALSE;
    }
    mass[0] = shape1->GetRandom();
    return kTRUE;
}

Double_t PMassSampling::GetWeight(Double_t *mass, Int_t *didx) {
    if (!shape1) {
	Warning("GetWeight","(%s): no TF1 found",GetDescription());
	return 0;
    }
    return shape1->Eval(mass[0]);
}

