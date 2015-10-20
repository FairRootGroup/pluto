// Author: I. Froehlich
// Written:2.6.2007
// Revised: 

#ifndef _PGENBOD_H_
#define _PGENBOD_H_



#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"

#define MAX_GENBOD_NUM 10


class PGenBod : public PChannelModel  {
  
 public:
    PGenBod();
    PGenBod(const Char_t *id, const Char_t *de, Int_t key);
    
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);
    Bool_t SampleMomentum(void);

    using PChannelModel::GetWeight;
    Double_t GetWeight();

    void UseWeights(void){SetVersionFlag(VERSION_WEIGHTING);};

    void   Print(const Option_t* delme=NULL) const;  //Debug info
    
 private:

    PParticle *parent, *primary, *daughter[MAX_GENBOD_NUM];
    int    n_daughters;
    Double_t local_weight,weight_max;

    PChannelModel *correlator;
    
    ClassDef(PGenBod,0) //Pluto's adaption from the original genbod
};

#endif


