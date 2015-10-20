// Author: I. Froehlich
// Written: 27.9.2006
// Revised: 

#ifndef _PSCATTERDISTRIBUTION_H_
#define _PSCATTERDISTRIBUTION_H_



#include "TF1.h"
#include "TF2.h"
#include "PDistribution.h"


class PScatterDistribution : public PDistribution  {
  
 public:
    PScatterDistribution();
    PScatterDistribution(const Char_t *id, const Char_t *de);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);

    Bool_t IsValid(void);


    void   SetAngleFunction(TF1 * f) {angles1=f;};
    void   SetAngleFunction(TF2 * f) {angles2=f;};


 private:
  
    TF1 * angles1;
    TF2 * angles2;
    PParticle * primary;
    PParticle * parent;
    PParticle * beam;
    PParticle * target;
    PParticle * mass_reference;

    ClassDef(PScatterDistribution,0)
};

#endif


