// Author: I. Froehlich
// Written: 16.9.2009
// Revised: 
// 

#ifndef _PETADOUBLEDALITZ_H_
#define _PETADOUBLEDALITZ_H_


#include "PChannelModel.h"
#include "TRandom2.h"




class PEtaDoubleDalitz : public PChannelModel  {
  
public:
  
  PEtaDoubleDalitz(const Char_t *id, const Char_t *de, Int_t key);
  PDistribution* Clone(const char*delme=NULL) const;
  
  Bool_t Init(void);
  
  using PChannelModel::SampleMass;
  Bool_t SampleMass(void);
  
private:
  
  Double_t   Gen2lepton1(Double_t m);
  TRandom2  *gRand;
  PParticle *dil1,*dil2,*parent;

  Double_t       ff_w_max;             //Max weight of FF model
  PChannelModel *formfactor_model;     //form factor object
  
  ClassDef(PEtaDoubleDalitz,0)  // Simple Eta Dalitz decay based on two dileptons

};

#endif // _PETADOUBLEDALITZ_H_
