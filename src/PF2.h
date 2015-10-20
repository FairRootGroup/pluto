// Author: I. Froehlich
// Written: 18.8.2010
// Revised: 

#ifndef _PF2_H_
#define _PF2_H_

#include "TF2.h"


class PF2: public TF2 {

 public:

    void SetEpsilon(Double_t e) {
	epsilon=e;
    };

    Bool_t MakeIntegral(void);
    PF2(const char *name, Double_t (*fcn)(Double_t *, Double_t *), 
	Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Int_t npar);

 protected:

    using TF2::Integral;
    Double_t Integral(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t epsilon);

    Double_t epsilon;

    ClassDef(PF2,0)  // Adapted TF2 class

};

#endif
