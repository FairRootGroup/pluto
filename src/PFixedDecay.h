// Author: I. Froehlich
// Written: 27.5.2007
// Revised: 

#ifndef _PFIXEDDECAY_H_
#define _PFIXEDDECAY_H_

#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"
#include "PDynamicData.h"
#include "PKinematics.h"

#define MAX_FIXED_NUM 10

class PFixedDecay : public PChannelModel  {
  
 public:
    PFixedDecay();
    PFixedDecay(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t   Init(void);
    Bool_t   SampleMass(void);
    Bool_t   SampleMass(Double_t *mass, Int_t *didx=NULL);

    int      GetDepth(int i=0);

    Bool_t   GetWidth(Double_t mass, Double_t *width, Int_t didx=-1);

 protected:

    PParticle *parent, *daughter[MAX_FIXED_NUM];
    int    parent_id,d_id[MAX_FIXED_NUM];          //PIDs
    double parent_mass, dmass[MAX_FIXED_NUM]; //Static Masses
    int    n_daughters;

    ClassDef(PFixedDecay,0)  // Decay of A -> a+b+c+... with fixed product masses
};

#endif


