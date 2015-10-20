// Author: I. Froehlich
// Written: 23.5.2007
// Revised: 

#ifndef _PHADRONMODEL_H_
#define _PHADRONMODEL_H_

#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"
#include "PDynamicData.h"

class PHadronModel : public PChannelModel  {
  
 public:
    PHadronModel();
    PHadronModel(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t GetWidth(Double_t mass, Double_t *width, Int_t didx=-1);


 protected:
    Double_t global_weight_scaling;
    ClassDef(PHadronModel,0)  // Base class for the shapes of resonances
};

#endif


