// Author: I. Froehlich
// Written: 27.5.2007
// Revised: 

#ifndef _PEEDIRECTDECAY_H_
#define _PEEDIRECTDECAY_H_

#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"
#include "PDynamicData.h"
#include "PKinematics.h"

class PEEDirectDecay : public PChannelModel  {
  
 public:
    PEEDirectDecay();
    PEEDirectDecay(const Char_t *id, const Char_t *de, Int_t key=-1);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t   Init(void);
    Bool_t   SampleMass(void);
    Bool_t   SampleMass(Double_t *mass, Int_t *didx=NULL);

    Bool_t   GetWidth(Double_t mass, Double_t *width, Int_t didx=-1);

    using PDistribution::GetWeight;   
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);
    int      GetDepth(int i=0);

    virtual Double_t	Eval(Double_t x, Double_t y = 0, Double_t z = 0, Double_t t = 0) const;
    //TF1 wrapper

    void SetPiCutoff(int i){use_pi_cutoff=i;};
    void SetHadronicPS(int i){use_hadronic_ps=i;};

 protected:
  
    PParticle *parent,*e1,*e2;
    
    int use_pi_cutoff,use_hadronic_ps;
    int parent_id;
    double cv,mlep;

    ClassDef(PEEDirectDecay,0)  // Direct decays of vector mesons -> dilepton/dimuon
};

#endif


