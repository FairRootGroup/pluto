/////////////////////////////////////////////////////////////////////
//
// Decay Model in 1 unstable and N stable hadrons
// The unstable hadron is which has the largest width
// All other masses are set to the fixed value
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "PHadronDecayM1N.h"


ClassImp(PHadronDecayM1N)

PHadronDecayM1N::PHadronDecayM1N()  {
} ;

PHadronDecayM1N::PHadronDecayM1N(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de,key) {
    if (is_channel<0)
	Warning("PHadronDecayM1N","The model (%s) should be bound to CHANNELS only",de);
  
    //Get particles
    Int_t tid[11];
    tid[0]=10; 
    makeStaticData()->GetDecayModeByKey(primary_key,tid); // retrieve current mode info

    //Parent ALWAYS important (also for the inherited classes)
    parent_id   = makeStaticData()->GetDecayParentByKey(primary_key);
    parent_mass = makeStaticData()->GetParticleMass(parent_id);

    parent = NULL;
    daughter_pos=0;
    unstable_pos=-1;
    unstable_width=0;
    for (int i=0; i<M1N_MAX_DAUGHTERS;i++) {
	daughters[i]= NULL;
    }
    for (int i=0; i<tid[0];i++) {
	daughter_masses[i]=makeStaticData()->GetParticleMass(tid[i+1]);
	daughter_id[i]=tid[i+1];
	daughter_didx[i]=-1;
	if (makeStaticData()->GetParticleTotalWidth(tid[i+1]) > unstable_width) {
	    //define M1 as the most unstable one
	    unstable_pos   = i;
	    unstable_id    = tid[i+1];
	    unstable_width = makeStaticData()->GetParticleTotalWidth(tid[i+1]);
	}
    }
    mesh  = NULL;
    event = new TGenPhaseSpace();
    
    maxmesh = 300;
    //We scan the complete phase space 
    old_parent_mass=0;

    mesh = new PMesh(maxmesh-2,"mesh");
} ;

PDistribution* PHadronDecayM1N::Clone(const char*delme) const {
    return new PHadronDecayM1N((const PHadronDecayM1N &)* this);
};

void PHadronDecayM1N::UpdateMesh(void) {
    
    Double_t mmin = PData::LMass(unstable_id);
    Double_t mmax = PData::UMass(unstable_id);
    Double_t dm=(mmax-mmin)/(maxmesh-3.);          // mass increment
    for (int i=0;i<maxmesh-2;++i) {
	Double_t mm=mmin+i*dm;                     // current invariant mass
	
	daughter_masses[unstable_pos] = mm;
	event->SetDecay(*(TLorentzVector*)parent,daughter_pos, daughter_masses);
	Double_t phase_space = 1./event->GetWtMax();

	//cout << mm << ":" << phase_space << endl;

	Double_t model_weight = unstable_model->GetWeight(mm,&unstable_didx);
	
	mesh->SetNode(i,phase_space*model_weight); 

    }
    mesh->SetMin(mmin);                 // store mass threshold
    mesh->SetMax(mmax);                 // store mass ceiling

}

Bool_t PHadronDecayM1N::Init(void) {
    //Init function called once for each PChannel

    unstable_particle = NULL;

    //looking for parent. This is mandatory
    parent = GetParticle("parent");
    if (!parent) {
	Warning("Init","Parent not found");
	return kFALSE;
    }

    daughter_pos=0; //clear stuff because the Attach function makes a clone
    
    //Now get all daughters
    for (int i=0; i<M1N_MAX_DAUGHTERS;i++) {
	daughters[i]= GetParticle("daughter");
	if (daughters[i]) {
	    if (daughters[i]->ID() == unstable_id) {
		unstable_particle = daughters[i];
	    }
	    daughter_pos++;
	}
    }

    if (!unstable_particle) {
	Warning("Init","Unstable particle not found");
	return kFALSE;
    }

    return kTRUE;
}

Bool_t PHadronDecayM1N::Prepare(void) {
    //Things which might change during the eventloop
    
    unstable_didx = unstable_particle->GetDecayModeIndex(1);
    return kTRUE;
}
 
int PHadronDecayM1N::GetDepth(int i) {
    //check if we have models
    //This also initializes the sub-models

    Int_t a1;
    unstable_model = makeDynamicData()->GetParticleModel(unstable_id);
    if (unstable_model) a1 = unstable_model->GetDepth(i);

    //makeStaticData()->SetDecayEmin(is_channel, all_masses);
    return a1; 
}

void PHadronDecayM1N::SubPrint(Int_t opt) const {
    //Print sub-models
    
    if (unstable_model) {cout << " ";cout << unstable_model->GetDescription();}
}


Double_t PHadronDecayM1N::EvalPar(const Double_t *x, const Double_t *params) {
    return Eval(x[0]);
}

Double_t PHadronDecayM1N::Eval(Double_t x, Double_t y , Double_t z , Double_t t ) const
{
 
    
    if (mesh)  {
	return mesh->GetLinearIP(x);
    }
    return 0.;
}

Bool_t PHadronDecayM1N::SampleMass(void) {
    //Mass-sampling wrapper

    parent_mass = parent->M();

    if (fabs( parent_mass - old_parent_mass) > M1N_PARENT_GRID) UpdateMesh();

    Double_t new_mass = this->GetRandom();
    
    unstable_particle->SetM(new_mass);

    return kTRUE;
};



ClassImp(PHadronDecayM1N)
