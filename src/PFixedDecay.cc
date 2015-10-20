/////////////////////////////////////////////////////////////////////
//
//
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PFixedDecay.h"

PFixedDecay::PFixedDecay()  {
} ;

PFixedDecay::PFixedDecay(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de,key) {

    if (is_channel<0)
	Warning("PFixedDecay","The model (%s) should be bound to CHANNELS only",de);
  
    //Get particles
    Int_t tid[11];
    tid[0]=10; 
    makeStaticData()->GetDecayModeByKey(primary_key,tid); // retrieve current mode info

    //Parent ALWAYS important (also for the inherited classes)
    parent_id   = makeStaticData()->GetDecayParentByKey(primary_key);
    parent_mass = makeStaticData()->GetParticleMass(parent_id);

    n_daughters = tid[0];
    if (n_daughters>MAX_FIXED_NUM) {
	Warning("PFixedDecay","n_daughters>MAX_FIXED_NUM");
    }

    parent=NULL;
    for (int i=0;i<n_daughters;i++) {
	dmass[i]=makeStaticData()->GetParticleMass(tid[i+1]);
	d_id[i]=tid[i+1];
	daughter[i]=NULL;
    }

    version_flag |= VERSION_MASS_SAMPLING;  //Only one mass sampling in the PChannel

} ;

PDistribution* PFixedDecay::Clone(const char*delme) const {
    return new PFixedDecay((const PFixedDecay &)* this);
};

Bool_t PFixedDecay::Init(void) {
    //Init function called once for each PChannel
    
    parent = GetParticle("parent");
    if (!parent) {
	Warning("Init","Parent not found");
	return kFALSE;
    }
 
    for (int i=0;i<n_daughters;i++) {
	daughter[i]=GetParticle(makeStaticData()->GetParticleName(d_id[i]));
	if (!daughter[i]) {
	    Warning("Init","daughters not found");
	    return kFALSE;
	}
    }
    
    return kTRUE;
}



int PFixedDecay::GetDepth(int i) {
    
    double mymin=0.;
    for (int i=1;i<=n_daughters;i++) {
	mymin+=dmass[i-1];
    }

    makeStaticData()->SetDecayEmin(is_channel, mymin);
    return 0; //stable products -> depth is 0
}

Bool_t PFixedDecay::SampleMass(void) {
    //Mass-sampling wrapper

    for (int i=0;i<n_daughters;i++) {
	daughter[i]->SetM(dmass[i]);
    }

    return kTRUE;
};

Bool_t PFixedDecay::SampleMass(Double_t *mass, Int_t *didx) {
    //Not much to do here...
    //Since we have 2 stable products, for completeness
    //we reset the masses to the nominal value

    for (int i=1;i<=n_daughters;i++) {
	mass[i]=dmass[i-1];
    }
    return kTRUE;
};




Bool_t PFixedDecay::GetWidth(Double_t mass, Double_t *width, Int_t didx) {
    //Returns the fixed width for stable parents
    //and the p.s. for the decay unstable -> stable stable
    //decays > 2 body are not supported yet

    if ((n_daughters>2) || (makeStaticData()->GetParticleTotalWidth(parent_id) <  
			    (*unstable_width )) ) {
	*width = makeStaticData()->GetDecayPartialWidth(is_channel);
	return kTRUE;
    } 

    double pCms=PKinematics::pcms(mass,dmass[0],dmass[1]);
    double pole_pCms=PKinematics::pcms(parent_mass,dmass[0],dmass[1]);

    *width = (pCms/pole_pCms) * makeStaticData()->GetDecayPartialWidth(is_channel);
    return kTRUE;

};

ClassImp(PFixedDecay)
