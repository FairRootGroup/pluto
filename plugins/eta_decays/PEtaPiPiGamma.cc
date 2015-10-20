/////////////////////////////////////////////////////////////////////
// 
// Decay eta -> pi+ pi- gamma
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

#include "PEtaPiPiGamma.h"


ClassImp(PEtaPiPiGamma)

PEtaPiPiGamma::PEtaPiPiGamma()  {
} ;

PEtaPiPiGamma::PEtaPiPiGamma(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de,key) {

    if (is_channel<0)
	Warning("PEtaPiPiGamma","This model should be bound to CHANNELS only");

    weight_max = 0.0004;

} ;

PDistribution* PEtaPiPiGamma::Clone(const char*delme) const {
    return new PEtaPiPiGamma((const PEtaPiPiGamma &)* this);
};

Bool_t PEtaPiPiGamma::Init(void) {
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
    gamma = GetParticle("g");
    if (!gamma) {
	Warning("Init","Gamma not found");
	return kFALSE;
    }
    
    m_formfactor_model=
	GetSecondaryModel("formfactor");


    //Fold the 2 max values
    if (m_formfactor_model) {
	Info("Init","Found FF model");
	Double_t ff_w_max = m_formfactor_model->GetWeightMax();
	if (ff_w_max<0) {
	    Warning("Init","Max value of the FF model not initialized");
	} else {
	    weight_max *= ff_w_max;
	}

    }

    return kTRUE;
}


Double_t PEtaPiPiGamma::GetWeight(void) {
    //Returns the weight of ref [L1]

    Double_t e_tilde = gamma->E();

    TLorentzVector twopi = ((*(TLorentzVector*) pip)
			    + (*(TLorentzVector*) pim));

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
    Double_t p_star_value = p_star.P();
    Double_t s_pipi = twopi.M(); //to be consistent I used the mass (not squared) - see SECONDARY_MODELS

    Double_t weight = (
	e_tilde*e_tilde * p_star_value*p_star_value
	* pow(sin(theta_pi),2));

    if (m_formfactor_model) {
	weight *= m_formfactor_model->GetWeight(s_pipi);
    }

    return weight;

}

Bool_t PEtaPiPiGamma::IsValid(void) {
    //Use rejection mode...

    if (GetVersionFlag() & VERSION_WEIGHTING) return kTRUE; 
    //...but not if weighting enabled.

    Double_t weight = GetWeight();
    
    if (weight>weight_max) {
	weight_max = weight*1.1;
	Warning("IsValid","Weight > max, new max is %lf",weight_max);
    }

    if ((weight/weight_max)>PUtils::sampleFlat()) return kTRUE; // sample now distribution
    
    return kFALSE;


}

