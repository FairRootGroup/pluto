// Author: M.A. Kagarlis
// Written:  02.02.00
// Revised:  19.06.03  R.H.
// Revised:  26.06.07  IF
// PSaid Class Header

// R. Arndt's p+p elastic scattering parametrization
#ifndef _PSAID_H_
#define _PSAID_H_

#include <TGlobal.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TArrayD.h>
#include "PAngularDistribution.h"

class PSaid : public PAngularDistribution {

 public:
    
    PSaid();
    PSaid(const Char_t *id, const Char_t *de);
    
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);
    Bool_t Prepare(void);
    Bool_t IsValid(void);
    
    
    double SamplePolarAngle(double);


    ~PSaid();

 private:

    double told;
 

    int dim;
    //TRandom3 *REngine;                   // private M-Twistor random number engine
    TArrayD y, aa;                       // private arrays for the sampling algorithm
    void defaults();

    double dsdw(double, double);   // diff. cross section (mb/sr) by sc. angle (deg), Tlab (GeV)
    double sample(double, double); // sampling algorithm (rejection method)

    ClassDef(PSaid,0) //Pluto SAID Class

};
#endif // _PSAID_H_
