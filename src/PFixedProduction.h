// Author: I. Froehlich
// Written: 22.10.2010

#ifndef _PFIXEDPRODUCTION_H_
#define _PFIXEDPRODUCTION_H_

#include "PChannelModel.h"
#include "PDynamicData.h"
#include "PKinematics.h"

class PFixedProduction : public PChannelModel  {
  
 public:
    PFixedProduction();
    PFixedProduction(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t   Init(void);
    Bool_t   SampleMass(void);
    Bool_t   SampleMass(Double_t *mass, Int_t *didx=NULL);

    Bool_t   SampleMomentum(void);

    int      GetDepth(int i=0);



 protected:

    PParticle *parent, *daughter;
    int    parent_id,   d_id;          //PIDs
    double parent_mass, dmass;         //Static Masses

    ClassDef(PFixedProduction,0)  // Production of a+b -> X
};

#endif


