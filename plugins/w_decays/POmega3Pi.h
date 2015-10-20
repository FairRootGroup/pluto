// Author: I. Froehlich & T. Scheib
// Written: 17.01.2010
// Revised: 

#ifndef _POMEGA3PI_H_
#define _POMEGA3PI_H_

#include "TF1.h"
#include "TF2.h"
#include "PChannelModel.h"
#include "PDynamicData.h"
#include "PKinematics.h"
#include "PPropagator.h"

class POmega3Pi : public PDistribution  {
  
 public:	

    POmega3Pi();
    POmega3Pi(const Char_t *id, const Char_t *de);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init(void); 
    Bool_t Prepare(void);
    Bool_t Finalize(void);
    Bool_t IsValid(void);
    Bool_t CheckAbort(void);

    void SetMax(Double_t x) {
	max = x;
    };

    
 private: 		

    PParticle * side_particle[2]; //2 additional particles
    PParticle * primary;
    PParticle * parent;

    Double_t max;
	
    double diffgam(double M00, double M01);
    PChannelModel* RhoPropagator;

    ClassDef(POmega3Pi,0)  //omega -> pi+pi-pi0

};


#endif


