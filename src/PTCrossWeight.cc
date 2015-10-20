/////////////////////////////////////////////////////////////////////
//
// Calculates the total cross section weights for 
// a decay (a+b) -> something
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PTCrossWeight.h"


ClassImp(PTCrossWeight)

PTCrossWeight::PTCrossWeight() {

} ;

PTCrossWeight::PTCrossWeight(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel (id, de,key) {

    tcross = NULL;
    tcrossg = NULL;
    tcrossc = 1.;
    spline=kFALSE;
    g_spline=NULL;
    scaling=1.;
    isExcessEnergy=isMevEnergy=kFALSE;

    daughter_pos=0;
    for (int i=0; i<TCROSS_MAX_DAUGHTERS;i++) {
	daughters[i]= NULL;
    }
    EnableWeighting();
    SetExpectedWeightMean(-1); //switch off re-normalization
} ;

PDistribution* PTCrossWeight::Clone(const char*delme) const {
    return new PTCrossWeight((const PTCrossWeight &)* this);
};

Bool_t PTCrossWeight::Init(void) {

    //Get beam and target
    beam = GetParticle("beam");
    target = GetParticle("target");

    parent = GetParticle("parent");
//     if (!beam || !target) {
// 	Warning("Init","Beam/target not found");
// 	return kFALSE;
//     }

    //Now get all daugthers
    for (int i=0; i<TCROSS_MAX_DAUGHTERS;i++) {
	daughters[i]= GetParticle("daughter");
	if (daughters[i]) {	    
	    daughter_pos++;
	}
    }
  
//    if (!tcross && !tcrossg) Fatal("Init","No function found");

    return kTRUE;    
};


Double_t PTCrossWeight::GetWeight(void) {

    Double_t q=parent->M();
    if (beam && target)
	q=( *(TLorentzVector*)beam + *(TLorentzVector*)target  ).M();


    if (isExcessEnergy) {
	for (int i=0; i<daughter_pos; i++)
	    q-= daughters[i]->M();
    }
    
    if (isMevEnergy) q*=1000.;
    Double_t ret=GetWeight(&q);


//    cout << "Q:"  << q << ":" << ret << endl;
    return ret;
    
};

Double_t PTCrossWeight:: GetWeight(Double_t *mass, Int_t * ) {
    Double_t ret=0;
    if (tcrossg) {
	ret=tcrossg->Eval(mass[0],g_spline,"");
	//We have to check that we are not below the first and the last point
	if (mass[0]<tcrossg->GetX()[0]) {
	    ret=tcrossg->GetY()[0];
	}
	if (mass[0]>((tcrossg->GetX())[tcrossg->GetN()-1])) {
	  ret=(tcrossg->GetY())[tcrossg->GetN()-1];
	}
	

    } else if (tcross)
	ret=tcross->Eval(mass[0]);
    else 
	ret=tcrossc;
    // cout << "ret:" << ret << endl;
    if (ret>0) return ret*scaling ;
    return 0.;

}
