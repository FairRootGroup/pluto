// Author: I. Froehlich
// Written: 27.5.2007
// Revised: 

#ifndef _PHADRONDECAYM3_H_
#define _PHADRONDECAYM3_H_

#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"
#include "PDynamicData.h"
#include "PKinematics.h"

class PHadronDecayM3 : public PChannelModel  {
  
 public:
    PHadronDecayM3();
    PHadronDecayM3(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;


    Bool_t Prepare(void);
    Bool_t Init(void);
    Bool_t SampleMass(void);
    Bool_t SampleMass(Double_t *mass, Int_t *didx=NULL);

    Bool_t   GetWidth(Double_t mass, Double_t *width, Int_t didx=-1);

    Double_t GetWeight(void);
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);
    int      GetDepth(int i=0);

    Double_t Eval(Double_t x, Double_t y = 0, Double_t z = 0, Double_t t = 0) const;
    //TF1 wrapper

    void SubPrint(Int_t opt) const ;  //Print sub-models

 protected:

    Double_t parent_mass, mass1, mass2, mass3;
    int parent_id, id1, id2, id3;
    PChannelModel *model1, *model2, *model3;
    PParticle *parent,*daughter1,*daughter2,*daughter3;
    int didx1,didx2,didx3;

    ClassDef(PHadronDecayM3,0)  // Decay of Hadron -> 3*Hadron
};

#endif


