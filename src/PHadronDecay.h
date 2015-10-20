// Author: I. Froehlich
// Written: 27.5.2007
// Revised: 

#ifndef _PHADRONDECAY_H_
#define _PHADRONDECAY_H_

#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"
#include "PDynamicData.h"
#include "PKinematics.h"

class PHadronDecay : public PChannelModel  {
  
 public:
    PHadronDecay();
    PHadronDecay(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t   Init(void);
    Bool_t   SampleMass(void);
    Bool_t   SampleMass(Double_t *mass, Int_t *didx=NULL);

    Bool_t   GetWidth(Double_t mass, Double_t *width, Int_t didx=-1);

    int      GetDepth(int i=0);
//    Bool_t GetBR(Double_t mass, Double_t *br, Double_t totalwidth=-1);

    Double_t Eval(Double_t x, Double_t y = 0, Double_t z = 0, Double_t t = 0) const;
    Double_t EvalPar(const Double_t *x, const Double_t *params);
    //TF1 wrapper

    void Use_m0_over_m(int i) {use_m0_over_m=i;};

    void SetCutoffVersion(int i) {cutoff_version=i;};

 protected:
    int    parent_id,id1,id2;          //PIDs
    double parent_mass, mass1, mass2;  //Static Masses
    double parent_g0; //Static Width of parent
    
    //Width Configuration
    int use_fixed_delta;
    double fixed_delta;
    int angular_l;
    int cutoff_l;
    int use_m0_over_m,cutoff_version; //0=theis,1=ernst,2=none
    double w0;

    double HadronWidth(const double & m, const double & ma, const double & mb);

    ClassDef(PHadronDecay,0)  // Decay of Hadron -> Hadron(stable) + Hadron(stable)
};

#endif


