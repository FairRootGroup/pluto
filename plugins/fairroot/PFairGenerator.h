// Author: Ingo Froehlich
// Written: 13/10/2009
// Modified: 
// PFairGenerator Class Header

#ifndef _PFAIRGENERATOR_H_
#define _PFAIRGENERATOR_H_

#define FAIRGENERATOR_STACKSIZE 1000

#include "PBulkInterface.h"

class PFairGenerator: public PBulkInterface {

 private:

    Int_t fNumberParticles;                           //number of stored pparticles
    PParticle *fLocalStack[FAIRGENERATOR_STACKSIZE];  //stack


 protected:
    

 public:

    PFairGenerator();

    bool Modify(PParticle ** stack, int *decay_done, int * num, int maxnum);  //get particle array

    Bool_t GetNextParticle(Int_t *pdgType, Double_t *px,  Double_t *py,  
			   Double_t *pz,  Double_t *vx,  Double_t *vy,  Double_t *vz);
    
    ClassDef(PFairGenerator,0) // Interface to FairRoots primary generator

};
#endif 

















