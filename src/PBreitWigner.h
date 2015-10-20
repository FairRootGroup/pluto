// Author: Kagarlis
// Written: 23.5.2007
// Revised: I. Froehlich (from PData)

#ifndef _PBREITWIGNER_H_
#define _PBREITWIGNER_H_



#include "TF1.h"
#include "TF2.h"
#include "PHadronModel.h"


class PBreitWigner : public PHadronModel  {
  
 public:
    PBreitWigner();
    PBreitWigner(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    using PDistribution::GetWeight;   
    virtual Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);
    using PDistribution::SampleMass;
    Bool_t SampleMass(Double_t *mass, Int_t *didx=NULL);

    void SetWidthModel(PChannelModel * m) {
	//overwrites the partial model from the data base
	width_model = m;
    };
    
 protected:
  
    Double_t mr; // Resonance pole mass
    PChannelModel * width_model;

    ClassDef(PBreitWigner,0)  // Breit Wigner with mass-dependent width
};

#endif


