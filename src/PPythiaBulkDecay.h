// Author: Ingo Froehlich
// Written: 10/07/2007
// Modified: 
// PPythiaBulkDecay Class Header

#ifndef _PPYTHIABULKDECAY_H_
#define _PPYTHIABULKDECAY_H_

class TPythia6;

#include "PBulkInterface.h"

class PPythiaBulkDecay: public PBulkInterface {

 private:

    double tauMax;       //  decay only particles with tau>tauMax
    TPythia6 *fPythia;   //! pointer to Pythia object
    

 protected:
    

 public:

  PPythiaBulkDecay();

  bool Modify(PParticle ** stack, int *decay_done, int * num, int maxnum);  //decay all particle
  
  void SetTauMax(double t){tauMax=t*1.e-9;}; // go to sec
  void SetPythia(TPythia6 *p) {fPythia=p;} // set pointer to Pythia
  
  ClassDef(PPythiaBulkDecay,0) // Let particles decay the Pythia way 
};
#endif 

















