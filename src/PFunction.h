// Author: I. Froehlich
// Written: 07.4.2008
// Revised: 

#ifndef _PFUNCTION_H_
#define _PFUNCTION_H_


#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "PChannelModel.h"

#include "PBatch.h"


class PFunction : public PChannelModel  {
  
 public:
    PFunction();
    PFunction(Char_t *id, Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);

    using PDistribution::GetWeight;   
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);	

    using TF1::SetFunction;
    void SetFunction(TF1 * my_tf1) {tf1=my_tf1;};
    void SetFunction(TF2 * my_tf2) {tf2=my_tf2;};
    void SetFunction(TF3 * my_tf3) {tf2=my_tf3;};
    void SetFunction(Double_t c) {constant=c;};
    Bool_t AddEquation(char * command);

    using PDistribution::SampleMass;   
    Bool_t SampleMass(Double_t *mass, Int_t *didx=NULL);

 private:

    TF1 * tf1;
    TF2 * tf2;
    TF3 * tf3;
    Double_t constant,*vf,*vx;
    PBatch *batch;

    ClassDef(PFunction,0)  //Function wrapper (e.g. for model pathes)

};

#endif


