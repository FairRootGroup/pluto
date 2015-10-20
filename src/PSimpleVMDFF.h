// Author: I. Froehlich
// Written: 09.11.2010
// Revised: 

#ifndef _PSIMPLEVMDFF_H_
#define _PSIMPLEVMDFF_H_

#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"
#include "PBatch.h"

class PSimpleVMDFF : public PChannelModel  {
  
 public:
    PSimpleVMDFF();
    PSimpleVMDFF(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);
    
    using PChannelModel::GetWeight;
	
    Double_t GetWeight(void);
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);

    void SetVectorMesonMass(Double_t x) {
	vector_meson_mass = x;
	vector_meson_mass2 = x*x;
    };

    Bool_t AddEquation(char * command);

 private:

    Double_t vector_meson_mass,vector_meson_mass2;  //value of m_V
    Double_t *vff2,*vq,*vq2;

    PBatch *batch;

    PParticle *dilepton, *dilepton2, *parent, *ep, *em;

    ClassDef(PSimpleVMDFF,0)  //Simple VMD FF

};

#endif


