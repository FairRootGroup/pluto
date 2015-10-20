/////////////////////////////////////////////////////////////////////
//
// Resonance mass sampling from a relativistic
// Breit-Wigner distribution with mass-dependent width Breit-Wigner
// (low-mass cutoff through Gamma0*(q/q0)**(2L+1) ).
//
//                                  Author:  Kagarlis 
//                                  Implemented: I. Froehlich
//                                  Written: 23.5.2007
//                                  Revised: 
// 
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PBreitWigner.h"


ClassImp(PBreitWigner)

PBreitWigner::PBreitWigner()  {

} ;

PBreitWigner::PBreitWigner(const Char_t *id, const Char_t *de, Int_t key) :
    PHadronModel(id, de, key) {

    width_model = NULL;
    mr=makeStaticData()->GetParticleMassByKey(key);

};

PDistribution* PBreitWigner::Clone(const char*delme) const {
    return new PBreitWigner((const PBreitWigner &)* this);
};

Bool_t PBreitWigner::SampleMass(Double_t *mass, Int_t *didx) {
    // Resonance mass sampling from a relativistic
    // Breit-Wigner distribution with mass-dependent width Breit-Wigner
    //(low-mass cutoff through Gamma0*(q/q0)**(2L+1) ).
    // Used e.g. for 3-body decays in PHadronDecayM3.

    if (didx) {
	if (didx_option != didx[0]) {
	    didx_option=didx[0];
	    SetParameter(1,didx_option);
	}
    } else {
	if (didx_option !=-1) {
	    didx_option=-1;
	    SetParameter(1,didx_option);
	}
    }
    

    mass[0] = this->GetRandom();

    return kTRUE;

}

Double_t PBreitWigner::GetWeight(Double_t *mass, Int_t *didx) {
// relativistic Breit-Wigner distribution function for particle "key"
    
    double m=mass[0];
    double m2=m*m, 	
	mm=mr*mr-m2;
    int didx_local=-1;
    if (didx) didx_local=didx[0];

    global_weight_scaling = makeDynamicData()->GetParticleScalingFactor(is_pid);

    double wmt = 1.;
    if (!GetWidth(m,&width)) {
	Warning("GetWeight","GetWidth failed");
	return -1;
    }
    //N.B.: If didx is used, we sample the partial decay width

    double g=width, 
	g2=g*g;
    double partial_width;
    if (didx_local>=0) {
	if (!GetWidth(m,&partial_width,didx_local)) {
	    Warning("GetWeight","GetWidth failed");
	    return -1;
	}	
	if (width_model) {
	    width_model->GetWidth(m,&partial_width,didx_local);
	}
    } else partial_width=width;


// I copied this comment from PData here for bookkeeping (IF)
//  if (tempMtScaling>0. && (id==41 || id==42 || id==43)) { // obsolete > v4.09
//    fold rho line shape with thermal distribution
//    double mcut = TMath::Max(m,0.28);  // cut mass at 2 pion masses
//    wmt = mtIntegral(mcut,tempMtScaling)/mtIntegral(Mass(id),tempMtScaling);
//  } 
//  return wmt*m2*g2/(mm*mm+m2*g2);  // ->  BW(Mres) = 1
//    cout << "result: " << global_weight_scaling*wmt*m2*partial_width/(mm*mm+m2*g2) << endl;

    Double_t w= global_weight_scaling*wmt*m2*partial_width/(mm*mm+m2*g2);

//    if (w==0 && (didx_local==96)) CRASH;
    return w;
}

