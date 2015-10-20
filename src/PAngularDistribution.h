// Author: I. Froehlich
// Written: 3.7.2006
// Revised: 

#ifndef _PANGULARDISTRIBUTION_H_
#define _PANGULARDISTRIBUTION_H_



#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TSpline.h"
#include "PDistribution.h"
#include "PKinematics.h"

#define MAX_ANG_NUM 10


class PAngularDistribution : public PDistribution  {
  
 public:
    PAngularDistribution();
    PAngularDistribution(const Char_t *id, const Char_t *de);
    
    PDistribution* Clone(const char*delme=NULL) const;
    
    Bool_t Init(void);
    Bool_t Prepare(void);
    Bool_t Finalize(void);
    Bool_t IsValid(void);
    Bool_t CheckAbort(void);
    Bool_t SampleAngle(void);

    void SetAngleFunction(TF1 * f) {angles1=f;};
    void SetAngleFunction(TF2 * f) {angles2=f;};
    void SetAngleFunction(TGraph * f, Bool_t useSpline=kFALSE, 
                          Bool_t useSymmetry=kFALSE) {
       anglesg=f;
       spline=useSpline;
       reflection_symmetry=useSymmetry;
       if (spline) g_spline=new TSpline3("",f);
    };

    void SetRotate(Bool_t t) {rotate=t;};
    void NeverAbort(Bool_t t) {never_abort=t;}; //Avoid to re-sample the complete reaction chain
    void ForceRejectionMethod(Bool_t t) {always_reject=t;};

    void Print(const Option_t* delme=NULL) const;  //Debug info

    virtual Double_t    Eval(Double_t x, Double_t y = 0, Double_t z = 0, Double_t t = 0) const;
    virtual Double_t    EvalPar(const Double_t *x, const Double_t *params);
    //TF1 wrapper

    virtual double SamplePolarAngle(double r=0);
    Bool_t direct_sampling_possible,direct_sampling_done; //Try to use direct sampling and not rejection

 protected:
  
    TF1 * angles1;
    TF2 * angles2;
    TGraph * anglesg;
    TSpline3 * g_spline;
    PParticle * reference;
    PParticle * base_reference;
    PParticle * primary;
    PParticle * parent;
    PParticle * align;
    PParticle * mass_reference, *daughter[MAX_ANG_NUM];
    PParticle * beam;
    PParticle * target;
    PParticle * ang_reference;
    int    n_daughters;
    Bool_t check_abort,never_abort,always_reject;
    Bool_t rotate;
    Bool_t align_is_daughter;
    Bool_t reflection_symmetry; //use symmetry around theta=90°
    Bool_t spline;

    PParticle primary_tmp, reference_tmp,
	ang_tmp;                          //for boosting and rotation
    Bool_t Rotate(Int_t);                 //Do temporary boosting and rotation
    Bool_t RotateBack(Int_t);             //Do temporary boosting and rotation

    Double_t q_value;

    Double_t weight_max;

    ClassDef(PAngularDistribution,0) //Multi-purpose angular distributions

};

#endif


