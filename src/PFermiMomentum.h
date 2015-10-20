// Author: I. Froehlich
// Written: 23.6.2007
// Revised: 

#ifndef _PFERMIMOMENTUM_H_
#define _PFERMIMOMENTUM_H_



#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TSpline.h"
#include "PChannelModel.h"
#include "PAdaptiveMeshN.h"


class PFermiMomentum : public PChannelModel  {
  
 public:
    PFermiMomentum();
    PFermiMomentum(const Char_t *id, const Char_t *de, Int_t key);
    ~PFermiMomentum();

    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);
    Bool_t Prepare(void);
    Bool_t IsValid(void);
    int GetDepth(int i);

    void   Print(const Option_t* delme=NULL) const;  //Debug info
    void   SubPrint(Int_t opt) const;
    
    using PDistribution::GetWeight;
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);

    Bool_t SampleMomentum(void);

    using PChannelModel::SampleMass;
    Bool_t SampleMass(void);

    void SetMomentumFunction(TGraph * f, Bool_t useSpline=kFALSE) {
       mom=f;
       spline=useSpline;
       if (spline) g_spline=new TSpline3("",f);
    };

    virtual Double_t	Eval(Double_t x, Double_t y = 0, Double_t z = 0, Double_t t = 0) const;
    virtual Double_t    EvalPar(const Double_t *x, const Double_t *params);
    //TF1 wrapper

 protected:
    
    double SampleFermi(double & px, double & py, double & pz);

 private:

    PParticle * beam;
    PParticle * target;
    PParticle * spectator;
    PParticle * parent;
    PParticle * p1,*p2;
    PParticle * composite;
    Int_t num_of_realevents,num_of_sampledevents;
    Int_t didx_composite;
    PChannelModel * composite_model;
    Int_t tcross_key;
    Int_t debug_print;
    Double_t px,py,pz,participant_mass;
    TLorentzVector my_beam;
    PAdaptiveMeshN * mesh;

    TGraph * mom;
    TSpline3 * g_spline;
    Bool_t spline;


    ClassDef(PFermiMomentum,0)  // Quasifree scattering using deuteron wave function
};

#endif


