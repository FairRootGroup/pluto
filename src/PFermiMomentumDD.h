// Author: I. Froehlich
// Written: 13.11.2007
// Revised: 

#include "PFermiMomentum.h"

//Class definition

class PFermiMomentumDD : public PFermiMomentum  {
  
 public:

    PFermiMomentumDD(Char_t *id, Char_t *de, Int_t key);
    
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);

    using PChannelModel::SampleMass;
    Bool_t SampleMass(void);

private:

    PParticle * beam;
    PParticle * target;
    PParticle * spectator1, *spectator2;
    PParticle * parent;
    PParticle * p1,*p2;
    PParticle * composite;

    ClassDef(PFermiMomentumDD,0) //N+N in d+d scattering

};



