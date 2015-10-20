// Author: Ingo Froehlich
// Written: 10/07/2007
// Modified: 
// PBulkInterface Class Header

#ifndef _PBULKINTERFACE_H_
#define _PBULKINTERFACE_H_

#include "PParticle.h"

#define PPROJECTOR_PRIORITY 99
#define FILTER_PRIORITY 60
#define DECAY_PRIORITY 50
#define FILEINPUT_PRIORITY 1


class PBulkInterface: public TObject {

 private:

    static Int_t gBulkCounter; 

    
 protected:
    
    Double_t current_weight;
    Int_t    bulk_id;      //Unique bulk ID
    Int_t    fPriority;    //order when adding the bulk in PReaction

 public:

    PBulkInterface();
    
    virtual bool Modify(PParticle ** array, int *decay_done, int * num, int maxnum);  //Modify particle bulk

    void SetWeight(Double_t c) {current_weight=c;};
    Double_t GetWeight(void) {return current_weight;};

    void SetPriority(Int_t p) {
	fPriority = p;
    };
    Int_t GetPriority() {return fPriority;};



ClassDef(PBulkInterface,0) // Pluto bulk interface base class
};


#endif 

















