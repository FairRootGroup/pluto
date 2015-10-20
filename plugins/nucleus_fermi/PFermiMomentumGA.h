// Author: L. Witthauer &
//         M. Dieterle
// Written 24.03.2009

#ifndef _PFERMIMOMENTUMGA_H_
#define _PFERMIMOMENTUMGA_H_

#include <iostream>
#include "PFermiDistributions.h"
#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"
#include "string.h"


class PFermiMomentumGA : public PChannelModel  {
   
 public:
    
    PFermiMomentumGA();
    PFermiMomentumGA(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char* delme=NULL) const;
    Bool_t Init(void);
    using PChannelModel::SampleMass;
    Bool_t SampleMass(void);
    Bool_t SampleMomentum(void);
    Double_t GetRandomFermiMomentum(Double_t & px, Double_t & py, Double_t & pz);
  
private:

    //    PFermiDistributions *Dist;   
    PParticle * beam;
    PParticle * target;
    PParticle * spectator;
    PParticle * participant;
    PParticle * parent;
    PParticle * p1,*p2;
    PParticle * composite;

    PChannelModel * fermi_model;


    ClassDef(PFermiMomentumGA,0) 

};

#endif
