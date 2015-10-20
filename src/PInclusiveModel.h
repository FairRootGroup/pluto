// Author: I. Froehlich
// Written: 27.9.2006
// Revised: 

#ifndef _PINCLUSIVEMODEL_H_
#define _PINCLUSIVEMODEL_H_

#define INCLUSIVE_MAX_DAUGHTERS 3

#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"

class PInclusiveModel : public PChannelModel  {
  
 public:

    PInclusiveModel();
    PInclusiveModel(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);
    using PChannelModel::SampleMass;

    using PChannelModel::GetWeight;
    Double_t GetWeight(void);
    
    Bool_t SampleMass(void);
    int    GetDepth(int i=0);

    void   SetSampleFunction(TF1 * f) {sample=f;};
    void   SetSampleFunction(Double_t f) {sample_d=f;};

 private:
    PParticle * parent, *daughters[INCLUSIVE_MAX_DAUGHTERS], *primary; 
    PChannelModel *daughter_models[INCLUSIVE_MAX_DAUGHTERS];
    int daughter_pos;
    TF1 * sample;
    Double_t sample_d;
    Double_t weight;

    void Setup();
  


    ClassDef(PInclusiveModel,0)  //Model for mass sampling for X in a->b+c+X
};



#endif
