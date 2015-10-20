/////////////////////////////////////////////////////////////////////
//
// Beam smearing models: Angular semaring is possible as well 
// as momentum smearing using TF1 objects. 
// Particles needed: beam and target (as grandparent), parent
// The particle from(beam, target) 
// with higher momentum is defined as the real beam
//
//                                  Author: I. Froehlich
//                                  Written: 18.7.2007
//                                  Revised: 
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PBeamSmearing.h"

ClassImp(PBeamSmearing)

PBeamSmearing::PBeamSmearing() {
} ;

PBeamSmearing::~PBeamSmearing() {
    if (mybeam) delete (mybeam);
    if (mytarget) delete (mytarget);
} ;


PBeamSmearing::PBeamSmearing(Char_t *id, Char_t *de) :
    PDistribution(id, de) {
    beam = NULL;
    mybeam = NULL;
    target = NULL;
    mytarget = NULL;
    parent = NULL;
    mom_smearing = NULL;
    ang_smearing = NULL;
    mom_function = NULL;
    thetaBeam=0.;
    phiBeam=0.;
    sigmaBeam=0.;
} ;

PDistribution* PBeamSmearing::Clone(const char*delme) const {
    return new PBeamSmearing((const PBeamSmearing &)* this);
};


Bool_t PBeamSmearing:: Init(void) {
    //Read beam, target and parent


    beam = GetParticle("beam");
    target = GetParticle("target");
    parent=GetParticle("parent");

    if (!beam || !target) {
	Error("Init","beam or target not found");
	return kFALSE;
    }

    if (!mybeam)
	mybeam = new PParticle(beam);
    else
	*mybeam = *beam;
    if (!mytarget)
	mytarget= new PParticle(target);
    else
	*mytarget = target;

    if (!parent){
	Error("Init","Parent not found");
	return kFALSE;
    }

    if (!mom_smearing && !ang_smearing) {
	Warning("Init","No smearing model found");
    }

    return kTRUE;
}

Bool_t PBeamSmearing:: Prepare(void) {
    //We use the prepare function for smearing since it might affect
    //the mass and momentum sampling done in the next steps

    Double_t angle=0;
    if (ang_smearing) angle=ang_smearing->GetRandom() *TMath::Pi() / 180.;

    Double_t mom=1.;
    if (mom_smearing) mom=mom_smearing->GetRandom();

    //restore saved particles
    *beam=*mybeam;
    *target=*mytarget;

    //Apply first the new function smearing
    //because the following step includes tilting
    if (beam->P() > target->P()) {
	Double_t mymom = beam->Rho();
	if (mom_function) mymom=mom_function->GetRandom();
	beam->SetMom(mymom*mom);
	beam->RotateX(angle);
	beam->RotateZ(PUtils::sampleFlat()*2*TMath::Pi());
    } else {
	Double_t mymom = target->Rho();
	if (mom_function) mymom=mom_function->GetRandom();
	target->SetMom(mymom*mom);
	target->RotateX(angle);
	target->RotateZ(PUtils::sampleFlat()*2*TMath::Pi());
    }

    //The following lines have been copied from the "old" PReaction.cc
    if (thetaBeam>0. || sigmaBeam>0.) {    // we have a skewed and/or smeared beam axis
	Double_t thB=0.;
	Double_t phB=0.;
	if (sigmaBeam>0.) {   // gaussian smearing of beam axis
	    //  thB = TMath::Abs(PUtils::sampleGaus(0.,sigmaBeam));
	    //  phB = 6.283185308*PUtils::sampleFlat();
	    Double_t thx = PUtils::sampleGaus(0.,sigmaBeam);   
	    // this gives better results
	    Double_t thy = PUtils::sampleGaus(0.,sigmaBeam);
	    thB = sqrt(thx*thx+thy*thy);
	    phB = TMath::ACos(thx/thB);
	    if (thy<0.) phB = 6.283185308-phB;
	}
	TVector3 beamAxis(TMath::Sin(thB)*TMath::Cos(phB),
			  TMath::Sin(thB)*TMath::Sin(phB),
			  TMath::Cos(thB));
	if (thetaBeam>0.) {   // skewing of beam axis (i.e. average beam is off z axis)
	    TVector3 skew(TMath::Sin(thetaBeam)*TMath::Cos(phiBeam),
			  TMath::Sin(thetaBeam)*TMath::Sin(phiBeam),
			  TMath::Cos(thetaBeam));
	    beamAxis.RotateUz(skew);
	}
	
	if (beam->P() > target->P()) {
	    beam->RotateUz(beamAxis);

	} else {
	    target->RotateUz(beamAxis);
	}

    }
    


    parent->Reconstruct();

//    parent->Print("scatter");

    return kTRUE;

}

void PBeamSmearing::SetReaction(char * reaction) {
    char *array[2];
    Int_t array_s=2; //max products
    PUtils::Tokenize(reaction, "+", array, &array_s);
    Add("q","parent");
    Add(array[0],"grandparent","beam");
    Add(array[1],"grandparent","target");
    NoDaughters();

}

void PBeamSmearing::Print(const Option_t* delme) const{

    PDistribution::Print();
}


