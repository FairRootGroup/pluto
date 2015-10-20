// Author: I. Froehlich
// Written: 16.9.2009
// Revised: 
// 

#ifndef _PETADOUBLEDALITZFF_H_
#define _PETADOUBLEDALITZFF_H_


#include "PChannelModel.h"

//Additional weight for the eta double Dalitz FF

class PEtaDoubleDalitzFF : public PChannelModel  {
  
public:
  
  using PChannelModel::GetWeight;
  PEtaDoubleDalitzFF(Char_t *id, Char_t *de, Int_t key);
  PDistribution* Clone(const char*delme=NULL) const;
  Double_t GetWeight();
  
  void SetLambda(Double_t  val) {Lambda=val;};

  Bool_t Init(void);
  
 private:

  PParticle *dil1,*dil2,*parent;   //decay particles
  Double_t  Lambda;                //pole factor

  ClassDef(PEtaDoubleDalitzFF,0)  //Eta DD FF
};

#endif
