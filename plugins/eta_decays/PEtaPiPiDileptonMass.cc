/////////////////////////////////////////////////////////////////////
// 
// Decay eta -> pi+ pi- dilepton
// Model to sample the dilepton mass
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

#include "PEtaPiPiDileptonMass.h"


ClassImp(PEtaPiPiDileptonMass)

PEtaPiPiDileptonMass::PEtaPiPiDileptonMass()  {
} ;

PEtaPiPiDileptonMass::PEtaPiPiDileptonMass(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de,key) {

    if (is_channel<0)
	Warning("PEtaPiPiDileptonMass","This model should be bound to CHANNELS only");

    m_pi    = makeStaticData()->GetParticleMass("pi+");
    mass_e  = makeStaticData()->GetParticleMass("e-"); 
    mass_ee = 2.*mass_e;
    //it is important to avoid energy conservation violation, otherwise
    //genbod will be much slower:
    SetDynamicRange(mass_ee,makeStaticData()->GetParticleMass("eta")-2*m_pi);
} ;

PDistribution* PEtaPiPiDileptonMass::Clone(const char*delme) const {
    return new PEtaPiPiDileptonMass((const PEtaPiPiDileptonMass &)* this);
};

Bool_t PEtaPiPiDileptonMass::Init(void) {
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
    ep = GetParticle("e+");
    if (!ep) {
	Warning("Init","e+ not found");
	return kFALSE;
    }
    em = GetParticle("e-");
    if (!em) {
	Warning("Init","e- not found");
	return kFALSE;
    }

    vmd_formfactor_model=
	GetSecondaryModel("formfactor");

    return kTRUE;
}


Double_t PEtaPiPiDileptonMass::GetWeight(Double_t *mass, Int_t * ) {

    //This is the mass-dependent part of eqn. 64 of Wirzba's report

    return GetMassWeight(mass[0]);
}

Double_t PEtaPiPiDileptonMass::GetMassWeight(Double_t mass) const {
    
    //Mass-dependent part of the eqn. 64

    if (mass < mass_ee) return 0.;

    Double_t weight = 1./(8*mass*mass);

    if (vmd_formfactor_model) {
	weight *= vmd_formfactor_model->GetWeight(mass);
    }

    return weight;

}

Bool_t PEtaPiPiDileptonMass::SampleMass(Double_t *mass, Int_t *didx) {

    mass[0]=this->GetRandom();
    return kTRUE;
}

Double_t PEtaPiPiDileptonMass::Eval(Double_t x, Double_t y , Double_t z , Double_t t ) const {

    //return mass of the dilepton for GetRandom sampling

    return GetMassWeight(x);

}

// Double_t PEtaPiPiDileptonMass::GetWeight(void) {

//     return GetWeight(dilepton->M());
// }

