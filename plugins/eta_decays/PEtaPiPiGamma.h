// Author: I. Froehlich
// Written: 01.10.2009
// Revised: 

#ifndef _PETAPIPIGAMMA_H_
#define _PETAPIPIGAMMA_H_

#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"
#include "PDynamicData.h"
#include "PKinematics.h"

class PEtaPiPiGamma : public PChannelModel  {
  
 public:
    PEtaPiPiGamma();
    PEtaPiPiGamma(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);

    using PDistribution::GetWeight;
    using PChannelModel::GetWeight;

    Double_t GetWeight(void);

    Bool_t IsValid(void);

 protected:
  
    PParticle *parent, *pip, *pim, *gamma;

    PChannelModel *m_formfactor_model;     //form factor object

    Double_t weight_max;               //Maximum
  
    ClassDef(PEtaPiPiGamma,0)  // Decay eta -> pi+ pi- gamma
};

#endif


