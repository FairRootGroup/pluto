/////////////////////////////////////////////////////////////////////
// 
// Decay eta -> pi+ pi- dilepton
// PiPi correlation including the angles after the dilepton sampling
//
// In total, the equation reads as follows:
//
//
//
// References:
// [L1] Thimo Petri and Andreas Wirzba
//      Internal Report, Juelich
// 
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PEtaPiPiDilepton.h"


ClassImp(PEtaPiPiDilepton)

PEtaPiPiDilepton::PEtaPiPiDilepton()  {
} ;

PEtaPiPiDilepton::PEtaPiPiDilepton(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de,key) {

    if (is_channel<0)
	Warning("PEtaPiPiDilepton","This model should be bound to CHANNELS only");

    m_pi    = makeStaticData()->GetParticleMass("pi+");
    mass_e  = makeStaticData()->GetParticleMass("e-"); 
    mass_eta= makeStaticData()->GetParticleMass("eta"); 
    mass_ee = 2.*mass_e;
    weight_max = 0.00000001;

} ;

PDistribution* PEtaPiPiDilepton::Clone(const char*delme) const {
    return new PEtaPiPiDilepton((const PEtaPiPiDilepton &)* this);
};

Bool_t PEtaPiPiDilepton::Init(void) {
    //Init function called once for each PChannel
    
    //looking for parent. This is mandatory
    parent = GetParticle("parent");
    if (!parent) {
	Warning("Init","Parent not found");
	return kFALSE;
    }
    pip = GetParticle("pi+");
    if (!pip) {
	Warning("Init","Pi+ not found");
	return kFALSE;
    }
    pim = GetParticle("pi-");
    if (!pim) {
	Warning("Init","Pi- not found");
	return kFALSE;
    }
    em = GetParticle("e-");
    if (!em) {
	Warning("Init","e- not found");
	return kFALSE;
    }
    ep = GetParticle("e+");
    if (!ep) {
	Warning("Init","e+ not found");
	return kFALSE;
    }
    return kTRUE;
}


Double_t PEtaPiPiDilepton::GetWeight(Double_t *mass, Int_t * ) {
    //Inputs: s_pipi, q, theta_pip, theta_e, phi
    //This is the second part of eqn. 64 of Wirzba's report

    Double_t s_pipi    = mass[0]*mass[0];
    Double_t q         = mass[1];
    Double_t q2        = q*q;
    Double_t theta_pip = mass[2];
    Double_t theta_ep  = mass[4];
    Double_t phi       = mass[5];

    Double_t sin_theta_pip   = sin(theta_pip);
    Double_t beta_pi2  = (1-((4*m_pi*m_pi)/s_pipi));
    Double_t weight = (s_pipi*beta_pi2 *sin_theta_pip*sin_theta_pip);

    //for the curly bracket only the first part is supported at the moment
    Double_t beta_e2  = (1-((4*mass_e*mass_e)/q2));
    Double_t weight_part2 = PKinematics::lambda(mass_eta*mass_eta,s_pipi,q2) *
	(1+beta_e2*sin(theta_ep)*sin(theta_ep) * sin(phi)*sin(phi)  );

    return weight * weight_part2;
}

Bool_t   PEtaPiPiDilepton::SampleMass(void) {
    return kTRUE; //keep masses
}

Double_t PEtaPiPiDilepton::GetWeight(void) {
    //Returns the weight of ref [L1]

    TLorentzVector twopi = ((*(TLorentzVector*) pip)
			    + (*(TLorentzVector*) pim));
    TLorentzVector dil   = ((*(TLorentzVector*) ep)
			    + (*(TLorentzVector*) em));
    

    PParticle p_star(pip);
    //rotate such that twopi points to z-axis:
    Double_t Phi = twopi.Phi();
    Double_t Theta = twopi.Theta();
    twopi.RotateZ(-Phi);
    twopi.RotateY(-Theta);
    p_star.RotateZ(-Phi);
    p_star.RotateY(-Theta);
    p_star.Boost(-twopi.BoostVector());

    Double_t theta_pi = TMath::Pi() - p_star.Theta();

    Double_t s_pipi = twopi.M(); 

    PParticle ep_tmp(pip);
    Phi   = dil.Phi();
    Theta = dil.Theta();
    dil.RotateZ(-Phi);
    dil.RotateY(-Theta);
    ep_tmp.RotateZ(-Phi);
    ep_tmp.RotateY(-Theta);
    ep_tmp.Boost(-dil.BoostVector());

    Double_t phi=p_star.Phi() - ep_tmp.Phi();

    return GetWeight(s_pipi, dil.M(), theta_pi, ep_tmp.Theta(), phi);

}

Bool_t PEtaPiPiDilepton::IsValid(void) {
    //Use rejection mode...

    if (GetVersionFlag(VERSION_WEIGHTING)) return kTRUE; 
    //...but not if weighting enabled.

    Double_t weight = GetWeight();
    
    if (weight>weight_max) {
	Warning("IsValid","Weight (%lf) > max (%lf)",weight,weight_max);
	weight_max = weight*1.1;
    }

    if ((weight/weight_max)>PUtils::sampleFlat()) return kTRUE; // sample now distribution
    
    return kFALSE;
    

    //return kFALSE;


}

