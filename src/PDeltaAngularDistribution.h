// Author: I. Froehlich
// Written: 23.6.2007
// Revised: 

#ifndef _PDELTAANGULARDISTRIBUTION_H_
#define _PDELTAANGULARDISTRIBUTION_H_



#include "TF1.h"
#include "TF2.h"
#include "PAngularDistribution.h"


class PDeltaAngularDistribution : public PAngularDistribution  {
  
 public:
    PDeltaAngularDistribution();
    PDeltaAngularDistribution(const Char_t *id, const Char_t *de);
    
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);
    Bool_t Prepare(void);
    Bool_t IsValid(void);
    
    void Print(const Option_t* delme=NULL) const;  //Debug info
    
    void SetTerm(int i){use_term=i;};

    double SamplePolarAngle(double);

 private:
    double ds_dt(double cos_th_cm);

    void getNN_DeltaN_param();
    // parameters for N+N->N+Delta
    double lambda2, prefac, mres;

    int anisotropy;  
    
    int use_term;


    ClassDef(PDeltaAngularDistribution,0) //Angular distributions is the NN->N Delta reaction
};

#endif


