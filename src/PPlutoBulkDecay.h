// Author: Ingo Froehlich
// Written: 10/07/2007
// Modified: 
// PPlutoBulkDecay Class Header

#ifndef _PPLUTOBULKDECAY_H_
#define _PPLUTOBULKDECAY_H_

#include "PBulkInterface.h"

class PPlutoBulkDecay: public PBulkInterface {

 private:

    double tauMax;  //decay only particles with tau>tauMax
    int didx;
    int stackchannel;
    int recursiveMode;

 protected:
    

 public:

  PPlutoBulkDecay();

  bool Modify(PParticle ** stack, int *decay_done, int * num, int maxnum);  //decay all particle
  
  void SetTauMax(double t){tauMax=t*1.e-9;}; // go to sec

  void SetRecursiveMode(int t) {
      //0: no recursive decay
      //1: recursive decay enabled
      recursiveMode=t;
  };
  
  ClassDef(PPlutoBulkDecay,0) // Let particles decay the Pluto way 
};
#endif 

















