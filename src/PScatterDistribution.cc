/////////////////////////////////////////////////////////////////////
//
// Scattering angular distribution
// This algorithm evaluates the angular distribution of
// the particle "primary" which must be a "DAUGHTER"
//
// It works similar to AngularDistrubution, but the reference vector
// is the scattering OBE between the beam and the target.
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PScatterDistribution.h"
#include "PKinematics.h"


ClassImp(PScatterDistribution)

PScatterDistribution::PScatterDistribution()  {
} ;

PScatterDistribution::PScatterDistribution(const Char_t *id, const Char_t *de) :
    PDistribution(id, de) {

    angles1 = NULL;
    angles2 = NULL;
    beam = NULL;
    target = NULL;
    parent = NULL;

} ;

PDistribution* PScatterDistribution::Clone(const char*delme) const {
    return new PScatterDistribution((const PScatterDistribution &)* this);
};


Bool_t PScatterDistribution::Init(void) {

    //is the needed function set?
    if ((angles1 == NULL ) && (angles2 == NULL )) {
	Warning("Init","Angular distribution not found");
	return kFALSE;
    }

    //looking for primary. This is mandatory
    primary = GetParticle("primary");
    if (!primary) {
	Warning("Init","Primary not found");
	return kFALSE;
    }
    //primary->Print();

    beam = GetParticle("beam");
    target = GetParticle("target");

    //now get the parent
    for (int i=0; i<position; i++) {
	if (particle_flag[i] == PARTICLE_LIST_PARENT)
	    parent=particle[i];
    }

    if (!parent) {
	Warning("Init","Parent not found");
	return kFALSE;
    }
    return kTRUE;    
};



Bool_t PScatterDistribution::IsValid(void) {
    

   
    PParticle tmp_primary(primary);  //particle under investigation. Make better a copy
    PParticle tmp_parent(parent);  //particle under investigation. Make better a copy

    tmp_primary.Boost(parent->BoostVector());  // go back to lab frame

    //Thats OK if we are in t-channel:
    PParticle tmp_beam(beam);
    PParticle tmp_target(target);

    double sam = 0.5, sum = 1.;
    double t,u,tu,i;
    t = 0.5;
    tu= 0;


    if (parent->GetValue(T_MATRIX, &t) && 
	parent->GetValue(U_MATRIX, &u) &&
	parent->GetValue(TU_MATRIX, &tu)&&
	parent->GetValue(CHANNEL_POS, &i) ) {
	
//	cout << t << ":" << u << ":" << tu<< ":" << i << endl;
	
	//Calculate which process is more likely
	sum = t+u+tu;

	
    } 
    
    sam=PUtils::sampleFlat()*sum;
    
    if (sam>(t+0.5*tu)) { //u-channel      
      tmp_beam = *target;
      tmp_target = *beam;
    }

    tmp_parent.Boost(-tmp_target.BoostVector()); //go into target system
    tmp_primary.Boost(-tmp_target.BoostVector()); //go into target system
    tmp_beam.Boost(-tmp_target.BoostVector());

    //now we have to make sure that the "beam" is pointing to z direction
    Double_t Phi = tmp_beam.Phi();
    Double_t Theta = tmp_beam.Theta();
    tmp_parent.RotateZ(-Phi);
    tmp_parent.RotateY(-Theta);
    tmp_primary.RotateZ(-Phi);
    tmp_primary.RotateY(-Theta);

    //rotate such that parent (e.g. the Resonance) points to z
    Phi = tmp_parent.Phi();
    Theta = tmp_parent.Theta();
    tmp_parent.RotateZ(-Phi);
    tmp_parent.RotateY(-Theta);
    tmp_primary.RotateZ(-Phi);
    tmp_primary.RotateY(-Theta);

    //now we have to boost the primary
    tmp_primary.Boost(-tmp_parent.BoostVector());

    Double_t tmp_c0=cos(tmp_primary.Theta());
    Double_t f;

    

    if (angles1)
	f=angles1->Eval(tmp_c0);
    else 
	f=angles2->Eval(tmp_c0,tmp_parent.M());
    
    //cout << tmp_c0 << ":" << f << endl;
    if (f>PUtils::sampleFlat()) return kTRUE; // sample now angular distribution
    

    return kFALSE;

};


