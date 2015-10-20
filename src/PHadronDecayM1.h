// Author: I. Froehlich
// Written: 27.5.2007
// Revised: 

#ifndef _PHADRONDECAYM1_H_
#define _PHADRONDECAYM1_H_

#include "TF1.h"
#include "TF2.h"
#include "PHadronDecay.h"
#include "PDynamicData.h"
#include "PKinematics.h"

class PHadronDecayM1 : public PHadronDecay  {
  
 public:
    PHadronDecayM1();
    PHadronDecayM1(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);
    Bool_t Prepare(void);
    void   SubPrint(Int_t opt) const;
    Bool_t SampleMass(void);
    Bool_t SampleMass(Double_t *mass, Int_t *didx=NULL);

    Bool_t   GetWidth(Double_t mass, Double_t *width, Int_t didx=-1);

    Double_t GetWeight(void);
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);
    int      GetDepth(int i=0);

    Double_t Eval(Double_t x, Double_t y = 0, Double_t z = 0, Double_t t = 0) const;
    Double_t EvalPar(const Double_t *x, const Double_t *params);
    //TF1 wrapper

 protected:
    int unstable_id, stable_id;
    int unstable_position,stable_position;
    Double_t unstable_mass, unstable_ml; //Unstable mass pole, threshold
    Double_t unstable_dynamic_mass;      //sampled unstable mass
    Double_t stable_mass;
    bool sampleM1(const double & ecm);
    Double_t scale;
    double BWWeight(const int & i1, const double & m, 
		    const double & m1, const double & m2, int didx_local1=-1, 
		    int i2=0, int didx_local2=-1);
    double maxBWWeight(const int &i1, const double &m, const double &lower,
		       const double &upper, double &max1, const double &m2, int didx_local1=-1, 
		       int i2=0, int didx_local2=-1);


    PParticle *parent, *daughter1,*daughter2;
    PChannelModel *model1, *model2;
    int didx1, didx2, didx_unstable;

    ClassDef(PHadronDecayM1,0)  // Hadron decay in 1 unstable and 1 stable product
};

#endif


