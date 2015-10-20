// Author: I. Froehlich
// Written: 09.11.2010
// Revised: 

#ifndef _PETAPIPIDILEPTON_H_
#define _PETAPIPIDILEPTON_H_

#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"
#include "PDynamicData.h"
#include "PKinematics.h"

class PEtaPiPiDilepton : public PChannelModel  {
  
 public:
    PEtaPiPiDilepton();
    PEtaPiPiDilepton(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);

    using PDistribution::GetWeight;
    using PChannelModel::GetWeight;

    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);
    Double_t GetWeight(void);

    using PDistribution::SampleMass;
    using PChannelModel::SampleMass;
    Bool_t   SampleMass(void);


    Bool_t IsValid(void);

 protected:
  
    PParticle *parent, *pip, *pim, *ep, *em;
    Double_t  m_pi, mass_ee, mass_e, mass_eta;
    Double_t weight_max;               //Maximum
  
    ClassDef(PEtaPiPiDilepton,0)  // pipi correlation in eta -> pi+ pi- dilepton
};

#endif


