/////////////////////////////////////////////////////////////////////
//
// A wrapper to use root-functions as Pluto-models
//
// Can be used for simple models
// or to use the PAdaptiveMeshN algorithm
//
// Works also in batch mode via AddEquation
// input: _x, output: _f
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PFunction.h"


ClassImp(PFunction)

PFunction::PFunction() {

} ;

PFunction::PFunction(Char_t *id, Char_t *de, Int_t key) :PChannelModel (id, de,key) {

    tf1 = NULL;
    tf2 = NULL;
    tf3 = NULL;
    constant=0;
    batch   = NULL;

} ;

PDistribution* PFunction::Clone(const char*delme) const {
    return new PFunction((const PFunction &)* this);
};

Bool_t PFunction::Init(void) {

    return kTRUE;    
};

Bool_t PFunction::AddEquation(char * command) {
    if (!batch) batch = new PBatch();
    
    vx = makeStaticData()->GetBatchValue("_x"); 
    vf = makeStaticData()->GetBatchValue("_f"); 

    return batch->AddCommand(command);
};

Double_t PFunction::GetWeight(Double_t *mass, Int_t *didx) {

    if (batch) {
	
	*vx = mass[0];
	batch->Execute();
	return *vf;

    }

    if (tf1) return tf1->Eval(mass[0]);
    if (tf2) return tf2->Eval(mass[0],mass[1]);
    if (tf2) return tf3->Eval(mass[0],mass[1],mass[2]);

    return constant;
    
};


Bool_t PFunction::SampleMass(Double_t *mass, Int_t *didx) {

    if (tf1) {
	mass[0]=tf1->GetRandom();
	return kTRUE;
    }

    return kFALSE;

}
