/////////////////////////////////////////////////////////////////////
//
// Production of a+b -> X (so only 1 daughter)
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PFixedProduction.h"


PFixedProduction::PFixedProduction()  {
} ;

PFixedProduction::PFixedProduction(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de,key) {

    if (is_channel<0)
	Warning("PFixedProduction","The model (%s) should be bound to CHANNELS only",de);
    
    //Get particles
    Int_t tid[11];
    tid[0]=10; 
    makeStaticData()->GetDecayModeByKey(primary_key,tid); // retrieve current mode info

    //Parent ALWAYS important (also for the inherited classes)
    parent_id   = makeStaticData()->GetDecayParentByKey(primary_key);
    parent_mass = makeStaticData()->GetParticleMass(parent_id);

    if (tid[0]>1) {
	Warning("PFixedProduction","Only for 'decay's to 1 daughter");
    }

    dmass=makeStaticData()->GetParticleMass(tid[1]);
    d_id=tid[1];

    parent=NULL;
    daughter=NULL;

    version_flag |= VERSION_MASS_SAMPLING;  //Only one mass sampling in the PChannel

} ;

PDistribution* PFixedProduction::Clone(const char*delme) const {
    return new PFixedProduction((const PFixedProduction &)* this);
};

Bool_t PFixedProduction::Init(void) {
    //Init function called once for each PChannel
    
    parent = GetParticle("parent");
    if (!parent) {
	Warning("Init","Parent not found");
	return kFALSE;
    }
 
    daughter=GetParticle(makeStaticData()->GetParticleName(d_id));
	if (!daughter) {
	  Warning("Init","daughter not found");
	  return kFALSE;
	}
    return kTRUE;
}



int PFixedProduction::GetDepth(int i) {
    return 0; 
}

Bool_t PFixedProduction::SampleMass(void) {
    //Mass-sampling wrapper
    //mass=parent mass
    
    //check if we are in the allowed mass region
    if ((parent->M() < PData::LMass(d_id)) ||
	(parent->M() > PData::UMass(d_id))) return kFALSE;

    daughter->SetM(parent->M());

    return kTRUE;
};

Bool_t PFixedProduction::SampleMass(Double_t *mass, Int_t *didx) {
    //Not much to do here...
    //mass=parent mass
    
    mass[0]=parent_mass;
    
    return kTRUE;
};

Bool_t PFixedProduction::SampleMomentum(void) {

    daughter->SetPxPyPzE(0,0,0,daughter->M());

    return kTRUE;

};

ClassImp(PFixedProduction)
