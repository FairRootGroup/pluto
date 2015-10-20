// Author: I. Froehlich 
// Written: 5.6.2007
// Revised: 

#ifndef _PMASSAMPLING_H_
#define _PMASSAMPLING_H_



#include "TF1.h"
#include "TF2.h"
#include "PHadronModel.h"


class PMassSampling : public PHadronModel  {
  
 public:
    PMassSampling();
    PMassSampling(Char_t *id, Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    using PDistribution::GetWeight;   
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);

    using PDistribution::SampleMass;
    Bool_t   SampleMass(Double_t *mass, Int_t *didx=NULL);

    void     SetSamplingFunction(TF1 * f) {shape1=f;};

 private:
  
    TF1 *shape1;  //mass sampling shape as a function of resonance mass (=x)

    ClassDef(PMassSampling,0)  // General purpose mass sampling model
};

#endif


