// Author: I. Froehlich
// Written: 09.11.2010
// Revised: 

#ifndef _PETAPIPIDILEPTONMASS_H_
#define _PETAPIPIDILEPTONMASS_H_

#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"
#include "PDynamicData.h"
#include "PKinematics.h"

class PEtaPiPiDileptonMass : public PChannelModel  {
  
 public:
    PEtaPiPiDileptonMass();
    PEtaPiPiDileptonMass(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);

    using PDistribution::GetWeight;
    using PDistribution::SampleMass;
    Bool_t SampleMass(Double_t *mass, Int_t *didx=NULL);
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);

    Double_t GetMassWeight(Double_t mass) const;

    Double_t Eval(Double_t x, Double_t y , Double_t z , Double_t t ) const;

 protected:
  
    PParticle *parent, *pip, *pim, *ep, *em;
    Double_t  m_pi, mass_ee, mass_e;

    PChannelModel *vmd_formfactor_model;     //form factor object
  
    ClassDef(PEtaPiPiDileptonMass,0)  // Decay eta -> pi+ pi- dilepton (mass sampling)
};

#endif


