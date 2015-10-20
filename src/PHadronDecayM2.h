// Author: I. Froehlich
// Written: 25.7.2007
// Revised: 

#ifndef _PHADRONDECAYM2_H_
#define _PHADRONDECAYM2_H_

#include "TF1.h"
#include "TF2.h"
#include "PHadronDecayM1.h"
#include "PDynamicData.h"
#include "PKinematics.h"

class PHadronDecayM2 : public PHadronDecayM1  {
  
 public:
    PHadronDecayM2();
    PHadronDecayM2(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);
    Bool_t Prepare(void);
    void   SubPrint(Int_t opt) const;
    Bool_t SampleMass(void);
    Bool_t SampleMass(Double_t *mass, Int_t *didx=NULL);

    
    Bool_t   GetWidth(Double_t mass, Double_t *width, Int_t didx=-1);

    using PHadronDecayM1::GetWeight;   
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);
    int      GetDepth(int i=0);

 
 protected:

    Double_t dynamic_mass1,dynamic_mass2; //return value from sampleM2

    bool sampleM2(const double & ecm);

    ClassDef(PHadronDecayM2,0)  // Hadron decay in 2 unstable products
};

#endif


