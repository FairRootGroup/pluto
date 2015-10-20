/////////////////////////////////////////////////////////////////////
// 
// Decay w -> pi+ pi- pi0
//
// References:
// [L1] Hadronic three-body decays of light vector mesons.
//      S. Leupold, (Frankfurt U.) , M.F.M. Lutz, (Darmstadt, GSI)
//      Jul 2008. 7pp.
//      Published in Eur.Phys.J.A39:205-212,2009.
//      e-Print: arXiv:0807.4686 [hep-ph] 
// 
// Authors:  I. Froehlich, T. Scheib and S. Leupold
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "POmega3Pi.h"


ClassImp(POmega3Pi)


POmega3Pi::POmega3Pi() {

    primary = NULL;
    parent = NULL;
 
} ;

POmega3Pi::POmega3Pi(const Char_t *id, const Char_t *de) :
    PDistribution(id, de) {

    primary = NULL;
    parent = NULL;
    RhoPropagator = NULL;
    max = -1;
} ;

PDistribution* POmega3Pi::Clone(const char*delme) const {
    return new POmega3Pi((const POmega3Pi &)* this);
};

Bool_t POmega3Pi::Init(void) {

    
    //looking for primary. This is mandatory
    primary = GetParticle("primary");
    if (!primary) {
	Warning("Init","Primary not found");
	return kFALSE;
    }

    //now get the parent
    for (int i=0; i<position; i++) {
	if (particle_flag[i] == PARTICLE_LIST_PARENT)
	    parent=particle[i];
    }
    if (!parent) {
	Warning("Init","Parent not found");
	return kFALSE;
    }

    int side = 0;

    for (int i=0; i<position; i++) {
	if ((particle_flag[i] & PARTICLE_LIST_DAUGHTER) 
	    && (particle[i] != primary)) {
	    if (side == 2) {
		Warning("Init","More then 2 side particle found");
		return kFALSE;
	    }
	    side_particle[side] = particle[i];
	    side++;
	}
    }

    if (side != 2) {
	Warning("Init","Less then 2 side particle found");
	return kFALSE;
    }

    return kTRUE;    
};

Bool_t POmega3Pi::Prepare(void) {
    return kTRUE;
};

Bool_t POmega3Pi::Finalize(void) {
    return kTRUE;
};

Bool_t POmega3Pi::CheckAbort(void) {
    return kFALSE;
};

Bool_t POmega3Pi::IsValid(void) {

    TLorentzVector M00 = (*(TLorentzVector *) primary) + (*(TLorentzVector *) side_particle[0]);
    TLorentzVector M01 = (*(TLorentzVector *) primary) + (*(TLorentzVector *) side_particle[1]); // Definiere meinen TLorentzvector für pi0/pi- auf x-Achse

    double myvariable = M01.M2(); 
    double myvariable2 = M00.M2(); 
	
    double factor = POmega3Pi::diffgam(myvariable, myvariable2);
    

    if (factor>max) {
	Warning("IsValid","Dalitz factor %f > max %f",factor,max);
	if (max < 0) {
	    Warning("IsValid","Dalitz factor max not set");
	}
	max = factor*1.02;
    } 

    if (factor<0  ) {
	return kFALSE;
    } 

    if (((factor/max))>PUtils::sampleFlat()) {   
	return kTRUE; 
    }

    return kFALSE;
};


double POmega3Pi::diffgam(double M00, double M01) {
    //Written by S. Leupold
    double pi = TMath::Pi();

    // physical parameters m_V, h_P, h_A, f, pion mass m_pi, omega mass m_omega
    double mv, hp, ha, f, mpi, mom;
    // third combination of invariant masses (squared)
    double M02;
    // matrix element C_omega -> 3pi
    TComplex cc;
    // parts of cc
    TComplex h1, h2, h3, addh;
    // phase space, returned variable
    double p, dg;
    
    TComplex a,b,c,d,e,ff;

    //Take Pluto build-in rho propagator
    if (RhoPropagator == NULL) {
	RhoPropagator = makeDynamicData()->GetParticleSecondaryModel("rho0","propagator");
	if (RhoPropagator == NULL) 
            Fatal ("diffgam","RhoPropagator not defined");
    }

 
    mv = 0.78;
    hp = 0.304;
    ha = 2.1;
    f = 0.09;
    mpi = 0.14;
    mom = 0.78;
    
    M02 = -M00 - M01 + TMath::Power(mom,2) + 3*TMath::Power(mpi,2);

  
    p = -1./3.*((TMath::Power(M00,2)*M01)/4. + (M00*TMath::Power(M01,2))/4. 
		- (M00*M01*TMath::Power(mom,2))/4. - 
		(3*M00*M01*TMath::Power(mpi,2))/4. + (TMath::Power(mom,4)*TMath::Power(mpi,2))/4. - 
		(TMath::Power(mom,2)*TMath::Power(mpi,4))/2. + TMath::Power(mpi,6)/4.);
    
    
    
    a = M00 + TMath::Power(mom,2);
    h1 = RhoPropagator->GetAmplitude(&M00)*a.Re();
    
    b = M01 + TMath::Power(mom,2);
    h2 = RhoPropagator->GetAmplitude(&M01)*b.Re();
    
    c = M02 + TMath::Power(mom,2);
    h3 = RhoPropagator->GetAmplitude(&M02)*c.Re();
    
    d = h1+h2;
    addh = d+h3;
    
    e = ha*hp*mv / (4.*TMath::Power(f,3)*mom);
    cc = addh * e;
    
    ff = p*cc*TComplex::Conjugate(cc);
    dg = ff/(8.*pi*pi*pi*32.*mom*mom*mom);
    
    return dg;
};
