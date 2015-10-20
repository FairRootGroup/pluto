// Author: I. Froehlich
// Written: 30.3.2010
// Revised: 
// 

#ifndef _PETADOUBLEDALITZENV_H_
#define _PETADOUBLEDALITZENV_H_


#include "PChannelModel.h"
#include "TRandom2.h"




class PEtaDoubleDalitzEnv : public PChannelModel  {
  
public:
  
  PEtaDoubleDalitzEnv(const Char_t *id, const Char_t *de, Int_t key);
  PDistribution* Clone(const char*delme=NULL) const;
  
  Bool_t Init(void);

  Bool_t Finalize(void);

  Bool_t EndOfChain(void);
  
  
private:
  
  PParticle *dil1,*dil2,*parent,*ep1,*ep2,*em1,*em2;
  
  ClassDef(PEtaDoubleDalitzEnv,0)  //Complete Eta Dalitz decay, enveloped over the decey chain

};

#endif // _PETADOUBLEDALITZENV_H_
