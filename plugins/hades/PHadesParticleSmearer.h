// Author: 
// Written: 12.5.2009
// Revised: 
// PHadesParticleSmearer Class Header

#ifndef _PHADESPARTICLESMEARER_H_
#define _PHADESPARTICLESMEARER_H_

#include "TObject.h"
#include "PBulkInterface.h"
#include "PUtils.h"

class PHadesParticleSmearer: public PBulkInterface {
 
 public:
    PHadesParticleSmearer();
    

    Bool_t Modify(PParticle ** stack, int *decay_done, int * num, int stacksize);   
    
    void SetResolutionFactor(Double_t x) {
	resolution_factor=x;
    };

 protected:

    static double detrk4[7][5] ;
    static double multscrk4[7][5] ;
    static double detrk3[7][5] ;
    static double multscrk3[7][5] ;
    static double detkick[7][5] ;
    TRandom3 rand;


    double GetResolution(PParticle *p, Int_t iSetup=4);
    void   Smear(PParticle *p, double gamma);
    double resolution_factor ;
 
    
    ClassDef(PHadesParticleSmearer,0)  // Momentum smearing of HADES particle tracks

};

#endif
