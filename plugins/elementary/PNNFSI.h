// Author: I. Froehlich
// Written: 15.8.2008
// Revised: 

#ifndef _PNNFSI_H_
#define _PNNFSI_H_

#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"

#define NNFSI_MAX_DAUGHTERS 10

class PNNFSI : public PChannelModel  {
  
 public:

    PNNFSI();
    PNNFSI(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);
    using PChannelModel::SampleMass;

    using PChannelModel::GetWeight;
    Double_t GetWeight(void);
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);
    
    void SetA0R0(Double_t a, Double_t r) {
	a0=a;r0=r;
    };


 private:

    void Setup(void);
    PParticle *n1,*n2, *daughters[NNFSI_MAX_DAUGHTERS]; 
    Int_t daughter_pos,triplet,is_np;
    Double_t r0,a0,weight;

    ClassDef(PNNFSI,0)  // NN FSI using the Jost function
};



#endif
