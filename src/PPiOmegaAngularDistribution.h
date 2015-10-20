// Author: I. Froehlich
// Written: 23.6.2007
// Revised: 

#ifndef _PPIOMEGAANGULARDISTRIBUTION_H_
#define _PPIOMEGAANGULARDISTRIBUTION_H_

#define PI_OMEGA_piNNw    1
#define PI_OMEGA_piPPpiw  2
#define PI_OMEGA_piPDw    3

#include "TF1.h"
#include "TF2.h"
#include "PAngularDistribution.h"


class PPiOmegaAngularDistribution : public PAngularDistribution  {
  
 public:
    PPiOmegaAngularDistribution();
    PPiOmegaAngularDistribution(const Char_t *id, const Char_t *de);
    
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);
    Bool_t Prepare(void);
    Bool_t IsValid(void);
    
    void   Print(const Option_t* delme=NULL) const;  //Debug info
    
    double SamplePolarAngle(double);

    void   SetVersion(int i) {version=i;};

 private:
  
    double e_cm;
    int version;
    void getPiN_wN_param();
    // parameters for pi+N->N+w
    double PiN_w_h, PiN_w_y[5], PiN_w_slope[4], PiN_w_area[4];

    ClassDef(PPiOmegaAngularDistribution,0)  //Angular distributions in pi+p -> omega + X
};

#endif


