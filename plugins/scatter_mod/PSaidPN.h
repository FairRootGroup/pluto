// PSaidPN Class Header


#ifndef _PSAIDPN_H_
#define _PSAIDPN_H_

#include <TGlobal.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TArrayD.h>
#include "PAngularDistribution.h"

class PSaidPN : public PAngularDistribution {

 public:
    
    PSaidPN();
    PSaidPN(const Char_t *id, const Char_t *de);
    
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void);
    Bool_t Prepare(void);
    Bool_t IsValid(void);
    
    
    double SamplePolarAngle(double);

 private:

    double told;
    int dim;
    TArrayD y, aa;                 // private arrays for the sampling algorithm

    double dsdw(double, double);   // diff. cross section (mb/sr) by sc. angle (deg), Tlab (GeV)
    double sample(double, double); // sampling algorithm (rejection method)

    ClassDef(PSaidPN,0) //Pluto SAID Class (for p+n)

};
#endif // _PSAIPND_H_
