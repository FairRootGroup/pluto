/////////////////////////////////////////////////////////////////////
//
// Mass sampling of the remaining mass "X" in inclusive reactions
// like Q->a+b+X
// The model is based on a TF1 function which should be valid
// between [0,1], which corresponds to the free energy
// (so 1 is the maximum always)
// Sampling is done for the "primary" particle
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PInclusiveModel.h"
#include "PDynamicData.h"


ClassImp(PInclusiveModel)

PInclusiveModel::PInclusiveModel()  {

    Setup();

} ;

PInclusiveModel::PInclusiveModel(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de, key) {

    Setup();

} ;

void PInclusiveModel::Setup(void) {
    sample = NULL;
    sample_d = 0.;
    parent = NULL;
    daughter_pos=0;
    for (int i=0; i<INCLUSIVE_MAX_DAUGHTERS;i++) {
	daughters[i]= NULL;
	daughter_models[i]= NULL;
    }
}

PDistribution* PInclusiveModel::Clone(const char*delme) const {
    return new PInclusiveModel((const PInclusiveModel &)* this);
};


Bool_t PInclusiveModel::Init(void) {

    daughter_pos=0; //clear stuff because the Attach function makes a clone

    //is the needed function set?
//     if (sample == NULL ) {
// 	Warning("Init","Mass sampling distribution not found");
// 	return kFALSE;
//     }

    //looking for primary. This is mandatory
    primary = GetParticle("primary");
    if (!primary) {
	Warning("Init","Primary not found");
	return kFALSE;
    }

    parent = GetParticle("parent");
    if (!parent) {
	Warning("Init","Parent not found");
	return kFALSE;
    }

    //Now get all daughters
    for (int i=0; i<INCLUSIVE_MAX_DAUGHTERS;i++) {
	daughters[i]= GetParticle("daughter");
	if (daughters[i]) {	    
	    daughter_pos++;
	}
    }
//    cout <<"daughter_pos"<< daughter_pos << endl;

    return kTRUE;    
};

int PInclusiveModel::GetDepth(int i) {
    //check if we have models
    //This also initializes the sub-models
    //returns 0 in case of stable daughters BUGBUG->check all others

    Int_t a1=-1; //Default is stable decay products
    
    for (int i=0; i<daughter_pos;i++) {
	daughter_models[i]=makeDynamicData()->GetParticleModel(daughters[i]->ID());
	if (daughter_models[i]) {
	    int b1=daughter_models[i]->GetDepth(i);
	    a1=TMath::Max(a1,b1);
	}
    }

    return a1; 
}

Bool_t PInclusiveModel::SampleMass(void) {
    //In the sampling function, we first call the sub-models to
    
    Double_t mass, total_mass=0;
    for (int i=0; i<daughter_pos;i++) {
	if (daughter_models[i]) {
	    daughter_models[i]->SampleMass(&mass);
	    daughters[i]->SetM(mass);
	}
	total_mass+=daughters[i]->M();
	    
    }

    Double_t q = parent->M();
    if ((q-total_mass)<0) {
	Warning("SampleMass","Not enough energy");
	return kFALSE;
    }

    Double_t frac = sample_d;
    weight = 1.;

    if (sample) {
	
	frac=sample->GetRandom();
	weight = sample->Eval(frac);
    } 

    primary->SetM(frac*(q-total_mass));
    //cout << total_mass << ":" << q-total_mass  << ":" << frac*(q-total_mass) << endl;
    dynamic_range=q-total_mass;
    

    return kTRUE;
};

Double_t PInclusiveModel::GetWeight(void) {
    return weight;
    //BUGBUG: It can happen that sampleMass was not called
}


ClassImp(PInclusiveModel)

