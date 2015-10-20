/////////////////////////////////////////////////////////////////////
//
// Generic distribution for 3body decay
// used so far for eta -> pi+pi-pi0
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PDalitzDistribution.h"
#include "PDynamicData.h"

ClassImp(PDalitzDistribution)

PDalitzDistribution::PDalitzDistribution() {
    Fatal("PDalitzDistribution","Wrong ctor called");
} ;

PDalitzDistribution::PDalitzDistribution(const Char_t *id, const Char_t *de) :
    PDistribution(id, de) {

    primary = NULL;
    parent  = NULL;
    batch   = NULL;
 
} ;

PDistribution* PDalitzDistribution::Clone(const char*delme) const {
    return new PDalitzDistribution((const PDalitzDistribution &)* this);
};

Bool_t PDalitzDistribution::AddEquation(char * command) {
    if (!batch) batch = new PBatch();

    vprimary   = makeDynamicData()->GetBatchParticle("_primary"); 
    vs1        = makeDynamicData()->GetBatchParticle("_s1"); 
    vs2        = makeDynamicData()->GetBatchParticle("_s2"); 
 
    vf         = makeStaticData()->GetBatchValue("_f"); 

    return batch->AddCommand(command);
};

Bool_t PDalitzDistribution::Init(void) {

    
    //looking for primary. This is mandatory
    primary = GetParticle("primary");
    if (!primary) {
	Warning("Init","Primary not found");
	return kFALSE;
    }

    //now get the parent
    for (int i=0; i<position; i++) {
	if (particle_flag[i] == PARTICLE_LIST_PARENT)
	    parent=particle[i];
    }
    if (!parent) {
	Warning("Init","Parent not found");
	return kFALSE;
    }

    int side = 0;

    if ((side_particle[0] = GetParticle("s1"))) side++;
    if ((side_particle[1] = GetParticle("s2"))) side++;

    if (!side_particle[0] && (side_particle[0] = GetParticle("daughter"))) side++;
    if (!side_particle[1] && (side_particle[1] = GetParticle("daughter"))) side++;

    if (side != 2) {
	Warning("Init","Less or more than 2 additional found");
	return kFALSE;
    }

    return kTRUE;    
};

Bool_t PDalitzDistribution::Prepare(void) {
    return kTRUE;
};

Bool_t PDalitzDistribution::Finalize(void) {
    return kTRUE;
};

Bool_t PDalitzDistribution::CheckAbort(void) {
    return kFALSE;
};

Bool_t PDalitzDistribution::IsValid(void) {
    
    // eta -> pi+ pi- pi0
    // see e.g. Ref. 8 

    double factor = 1.;

    if (batch) {	

	*vprimary = primary;
	*vs1 = *side_particle[0];
	*vs2 = *side_particle[1];

	batch->Execute();
	factor =  *vf;
    } else {


	double eta_Q = parent->M() // eta
	    - primary->M()
	    - side_particle[0]->M()
	    - side_particle[1]->M();
	
	
	double pi0_T = primary->E() - primary->M();
	
	double dalitz_y = (3*pi0_T / eta_Q) -1;
	factor = 1 + slope1*dalitz_y + slope2*dalitz_y*dalitz_y;

    }
    
    if (factor>max) Warning("IsValid","Dalitz factor > max");
 
    if ((factor/max)>PUtils::sampleFlat()) return kTRUE; // sample now angular distribution
    
    return kFALSE;

};


