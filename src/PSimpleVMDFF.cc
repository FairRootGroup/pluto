/////////////////////////////////////////////////////////////////////
//
// Input:  the mass q (not squared - see SECONDARY_MODELS)
// Output: the form factor FF squared
//
// This class can be used in 2 modes:
//
// i)  either one has to set the vector meson mass,
//     the the following equation is used:
//       FF = m_V^2 / (m_V^2 - q^2)
//
// ii) or one can use a self-defined equation using the 
//     PBatch syntax
//     in this case the batch variables "_q" and "_q2"
//     are the dilepton mass (resp. squared)
//     The equation has to calculate _ff2, the form factor
//     squared. Example:
//
//     AddEquation("_ff2 = 0.17918225/( (0.4225- _q2)*(0.4225- _q2) + 0.000676)");
//
//     reproduces the Pluto build-in VMD function of Landsberg.
//
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PSimpleVMDFF.h"


ClassImp(PSimpleVMDFF)

PSimpleVMDFF::PSimpleVMDFF() {

} ;

PSimpleVMDFF::PSimpleVMDFF(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel (id, de,key) {

    dilepton = dilepton2 = parent = NULL;

    vector_meson_mass2 = -1;
    batch=NULL;

} ;

PDistribution* PSimpleVMDFF::Clone(const char*delme) const {
    return new PSimpleVMDFF((const PSimpleVMDFF &)* this);
};

Bool_t PSimpleVMDFF::AddEquation(char * command) {
    if (!batch) batch = new PBatch();

    vq   = makeStaticData()->GetBatchValue("_q"); 
    vq2  = makeStaticData()->GetBatchValue("_q2"); 
    vff2 = makeStaticData()->GetBatchValue("_ff2"); 

    return batch->AddCommand(command);

}

Bool_t PSimpleVMDFF::Init(void) {

    dilepton = GetParticle("dilepton");
    if (!dilepton) {
	dilepton = GetParticle("dimuon");
    }
    if (!dilepton) {
	ep = GetParticle("e+");
	em = GetParticle("e-");
    }
    if (!dilepton && !(ep && em) ) {
	Error("Init","No dilepton found");
	return kFALSE;
    }

    //double Dalitz?
    dilepton2 = GetParticle("dilepton");
    if (!dilepton2) {
	dilepton2 = GetParticle("dimuon");
    }

    parent = GetParticle("parent");
    if (!parent)  {
	    Error("Init","No parent found");
	    return kFALSE;
    }

    if (vector_meson_mass2<0 && !batch) {
	Error("Init","Vector meson mass OR eqn. must be initialized");
	return kFALSE;
    } 

    return kTRUE;    
};


Double_t PSimpleVMDFF::GetWeight(void) {

    Double_t q     = 0;
    if (dilepton)
	q     = dilepton->M();
    else {
	TLorentzVector dil = (*(TLorentzVector*)ep)+(*(TLorentzVector*)em);
	q     = dil.M();
    }
	
    Double_t pmass = parent->M();



    Double_t ff = GetWeight(q,pmass);

    if (dilepton2) {
	Double_t q     = dilepton2->M();
	ff *= GetWeight(q,pmass);
    }

    return ff;
    
};

Double_t PSimpleVMDFF:: GetWeight(Double_t *mass, Int_t * ) {
    Double_t q2 = mass[0]*mass[0];

    if (batch) {	
	*vq=mass[0];
	*vq2=q2;
	batch->Execute();
	return *vff2;
    }

    Double_t ff= vector_meson_mass2/(vector_meson_mass2 - q2);
    return ff*ff;

}
