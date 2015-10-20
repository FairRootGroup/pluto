// Author: Ingo Froehlich
// Written: 11/07/2007
// Modified: 
// PEmbeddedParticles Class Header

#ifndef _PEMBEDDEDPARTICLES_H_
#define _PEMBEDDEDPARTICLES_H_

#include "PBulkInterface.h"

#define MAX_EMBEDDEDPARTICLES 500

class PEmbeddedParticles: public PBulkInterface {

 private:

    PParticle *local      [MAX_EMBEDDEDPARTICLES];
    Int_t downscaling     [MAX_EMBEDDEDPARTICLES];
    Int_t last_downscaling[MAX_EMBEDDEDPARTICLES];
    Int_t local_pos;

    Int_t local_version[MAX_EMBEDDEDPARTICLES];

    //Cone sampling:
    Double_t local_pMin        [MAX_EMBEDDEDPARTICLES];
    Double_t local_pMax        [MAX_EMBEDDEDPARTICLES];
    Double_t local_mMin        [MAX_EMBEDDEDPARTICLES];
    Double_t local_mMax        [MAX_EMBEDDEDPARTICLES];
    Double_t local_openingAngle[MAX_EMBEDDEDPARTICLES];
    Double_t local_theta       [MAX_EMBEDDEDPARTICLES];
    Double_t local_phi         [MAX_EMBEDDEDPARTICLES];
    
    //Alternative "sector" sampling
    Double_t local_thetaMin[MAX_EMBEDDEDPARTICLES];
    Double_t local_thetaMax[MAX_EMBEDDEDPARTICLES];
    Double_t local_phiMin  [MAX_EMBEDDEDPARTICLES];
    Double_t local_phiMax  [MAX_EMBEDDEDPARTICLES];
    Int_t    nParticle;
    Double_t startPhi;

    //From batch system
    Double_t *vertex_x,*vertex_y,*vertex_z;

 protected:
    
    
 public:
    
    PEmbeddedParticles();
    
    Bool_t Modify(PParticle ** stack, int *decay_done, int * num, int stacksize);  //bulk interface
    
    Bool_t AddParticle(PParticle * particle, int downsc = 1);
    
    
    Bool_t SetSampling(Double_t pMin, Double_t pMax,
		       Double_t openingAngle, Double_t theta, Double_t phi,
		       Double_t mMin = 0., Double_t mMax = -.1);
    
    Bool_t SetSamplingSector(Double_t pMin    , Double_t pMax,
			     Double_t thetaMin, Double_t thetaMax,
			     Double_t phiMin  , Double_t phiMax,
			     Double_t phiStartVal, Int_t numParticle,
			     Double_t deltaPhi = 60., Int_t numSectors = 6);

    ClassDef(PEmbeddedParticles,0) // Add embedded particles in a PReaction
};
#endif 

















