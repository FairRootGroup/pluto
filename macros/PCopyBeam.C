#include "../src/PBulkInterface.h"

class PCopyBeam : public PBulkInterface {

 public:
    
    Bool_t Modify(PParticle ** stack, int *decay_done, int * num, int stacksize);  //bulk interface

    ClassDef(PCopyBeam,0) 
};



bool PCopyBeam::Modify(PParticle ** stack, int *decay_done, int * num, int stacksize) {

    for (int i=0; i< *num; i++) {
	PParticle * cur = stack[i];
	if ((cur->ID() > 1000) && (*num<(stacksize-2))) {
	    stack[*num] = cur->GetScattering(0);
	    stack[*num+1] = cur->GetScattering(1);
	    (*num) += 2;
	}
    }
    return kTRUE;
};

















