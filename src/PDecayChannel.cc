////////////////////////////////////////////////////////////////////
//  PDecayChannel
//
//  Author:   Volker Hejny
//  Written:  27.08.99
//  Revised:  30.06.00 MK
//  Revised:  23.07.07 IF (adaption for new framework)
//
//  PDecayChannel implements a linked list of decay modes. 
// 
////////////////////////////////////////////////////////////////////
using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>
#include "PDecayChannel.h"
#include "PParticle.h"

ClassImp(PDecayChannel)

// --------------------------------------------------------------------------
    PDecayChannel::PDecayChannel() {
    // 
    // Contructs an 'empty' channel (indicated by 'Daughters' containing
    // a NULL pointer).
    //  
    Weight = 0;
    NumOfDaughters = 0;
    Next = NULL;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
PDecayChannel::~PDecayChannel() {
    delete Next;
    Daughters.~TArrayI();
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
PDecayChannel* PDecayChannel::Unlink() {
    // 
    // The 'Next' element is set to NULL. The old value (i.e. the pointer
    // to the next possible decay mode) is returned.
    //
    PDecayChannel* temp = Next;
    Next = NULL;
    return temp;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel::AddChannel(PDecayChannel * newchannel) {
    // If the current channel is empty, it is filled with the contents of 
    // 'newchannel'. The object pointed to by 'newchannel' is deleted and
    // the value of 'newchannel' is replaced by a pointer to the current
    // channel.
    // If the current channel is not empty, but has no successor in 'Next',
    // 'Next' will point to 'newchannel'.
    // If the current channel contains a 'Next' object, the add request
    // is passed to this one. 
    // The other AddChannel functions work in the same way. The channel is
    // filled in the order current channel, 'Next' pointer and 'Next' object.
    if (!Daughters.GetArray()) {
	Double_t	w;
	Int_t       nd;
	Int_t*      ld;
	w  = newchannel->GetWeight();
	ld = newchannel->GetDaughters(nd);
	AddChannel(w,nd,ld);
	if (newchannel->GetNext()) AddChannel(newchannel->Unlink());
	delete newchannel;
	newchannel = this;
    }
    else {
	if (Next) Next->AddChannel(newchannel);
	else      Next = newchannel;
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel:: AddChannel(Double_t w, Int_t d1) {
    // The channel describes a decay into 1 daughter particle with particle id
    // d1 and branching ratio w. 
    // (This strange case is implemented only for completeness.)
    if (!Daughters.GetArray()) {
	Weight = w;
	NumOfDaughters = 1;
	Daughters.Set(1);
	Daughters[0] = d1;
    }
    else {
	if (!Next) Next = new PDecayChannel();
	Next->AddChannel(w,d1);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel:: AddChannel(Double_t w, Int_t d1, Int_t d2) {
    // The channel describes a decay into 2 daughter particles with particle 
    // ids d1, d2 and branching ratio w. 
    if (!Daughters.GetArray()) {
	Weight = w;
	NumOfDaughters = 2;
	Daughters.Set(2);
	Daughters[0] = d1;
	Daughters[1] = d2;
    }
    else {
	if (!Next) Next = new PDecayChannel();
	Next->AddChannel(w,d1,d2);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel:: AddChannel(Double_t w, Int_t d1, Int_t d2, Int_t d3) {
    // The channel describes a decay into 3 daughter particles with particle 
    // ids d1, d2, d3 and branching ratio w. 
    if (!Daughters.GetArray()) {
	Weight = w;
	NumOfDaughters = 3;
	Daughters.Set(3);
	Daughters[0] = d1;
	Daughters[1] = d2;
	Daughters[2] = d3;
    }
    else {
	if (!Next) Next = new PDecayChannel();
	Next->AddChannel(w,d1,d2,d3);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel:: AddChannel(Double_t w, Int_t d1, Int_t d2, Int_t d3, Int_t d4) {
    // The channel describes a decay into 4 daughter particles with particle 
    // ids d1, d2, d3, d4 and branching ratio w. 
    if (!Daughters.GetArray()) {
	Weight = w;
	NumOfDaughters = 4;
	Daughters.Set(4);
	Daughters[0] = d1;
	Daughters[1] = d2;
	Daughters[2] = d3;
	Daughters[3] = d4;
    }
    else {
	if (!Next) Next = new PDecayChannel();
	Next->AddChannel(w,d1,d2,d3,d4);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel:: AddChannel(Double_t w, Int_t nd, Int_t* ld) {  
    // The channel describes a decay into nd daughter particles with particle 
    // ids stored in the integer array ld and branching ratio w. 
    if (!Daughters.GetArray()) {
	Weight = w;
	NumOfDaughters = nd;
	Daughters.Set(nd);
	for (Int_t i=0; i<nd; i++) {
	    Daughters[i] = ld[i];
	}
    }
    else {
	if (!Next) Next = new PDecayChannel();
	Next->AddChannel(w,nd,ld);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel:: AddChannel(Double_t w, PParticle* d1) {
    // The channel describes a decay into 1 daughter particle represented by
    // a pointer to PParticle d1 and branching ratio w. 
    // (This strange case is implemented only for completeness.)
    if (!Daughters.GetArray()) {
	Weight = w;
	NumOfDaughters = 1;
	Daughters.Set(1);
	Daughters[0] = d1->ID();
    }
    else {
	if (!Next) Next = new PDecayChannel();
	Next->AddChannel(w,d1);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel:: AddChannel(Double_t w, PParticle* d1, PParticle* d2) {
    // The channel describes a decay into 2 daughter particles represented by
    // pointers to PParticle d1, d2 and branching ratio w. 
    if (!Daughters.GetArray()) {
	Weight = w;
	NumOfDaughters = 2;
	Daughters.Set(2);
	Daughters[0] = d1->ID();
	Daughters[1] = d2->ID();
    }
    else {
	if (!Next) Next = new PDecayChannel();
	Next->AddChannel(w,d1,d2);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel:: AddChannel(Double_t w, PParticle* d1, PParticle* d2, PParticle* d3) {
    // The channel describes a decay into 3 daughter particles represented by
    // pointers to PParticle d1, d2, d3 and branching ratio w. 
    if (!Daughters.GetArray()) {
	Weight = w;
	NumOfDaughters = 3;
	Daughters.Set(3);
	Daughters[0] = d1->ID();
	Daughters[1] = d2->ID();
	Daughters[2] = d3->ID();
    }
    else {
	if (!Next) Next = new PDecayChannel();
	Next->AddChannel(w,d1,d2,d3);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel:: AddChannel(Double_t w, PParticle* d1, PParticle* d2, PParticle* d3, PParticle* d4) {
    // The channel describes a decay into 4 daughter particles represented by
    // pointers to PParticle d1, d2, d3, d4 and branching ratio w. 
    if (!Daughters.GetArray()) {
	Weight = w;
	NumOfDaughters = 4;
	Daughters.Set(4);
	Daughters[0] = d1->ID();
	Daughters[1] = d2->ID();
	Daughters[2] = d3->ID();
	Daughters[3] = d4->ID();
    }
    else {
	if (!Next) Next = new PDecayChannel();
	Next->AddChannel(w,d1,d2,d3,d4);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel:: AddChannel(Double_t w, Int_t nd, PParticle** ld) {  
    // The channel describes a decay into nd daughter particles represented by
    // pointers to PParticle stored in the array ld and branching ratio w. 
    if (!Daughters.GetArray()) {
	Weight = w;
	NumOfDaughters = nd;
	Daughters.Set(nd);
	for (Int_t i=0; i<nd; i++) {
	    Daughters[i] = ld[i]->ID();
	}
    }
    else {
	if (!Next) Next = new PDecayChannel();
	Next->AddChannel(w,nd,ld);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel:: AddChannel(Double_t w, char* d1) {
    // The channel describes a decay into 1 daughter particle represented by
    // its name d1 and branching ratio w.
    // (This strange case is implemented only for completeness.)
    if (!Daughters.GetArray()) {
	Weight = w;
	NumOfDaughters = 1;
	Daughters.Set(1);
	Daughters[0] =makeStaticData()->GetParticleID(d1);
    }
    else {
	if (!Next) Next = new PDecayChannel();
	Next->AddChannel(w,d1);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel:: AddChannel(Double_t w, char* d1, char* d2) {
    // The channel describes a decay into 2 daughter particles represented by
    // their names d1, d2 and branching ratio w.
    if (!Daughters.GetArray()) {
	Weight = w;
	NumOfDaughters = 2;
	Daughters.Set(2);
	Daughters[0] =makeStaticData()->GetParticleID(d1);
	Daughters[1] =makeStaticData()->GetParticleID(d2);
    }
    else {
	if (!Next) Next = new PDecayChannel();
	Next->AddChannel(w,d1,d2);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel:: AddChannel(Double_t w, char* d1, char* d2, char* d3) {
    // The channel describes a decay into 3 daughter particles represented by
    // their names d1, d2, d3 and branching ratio w.
    if (!Daughters.GetArray()) {
	Weight = w;
	NumOfDaughters = 3;
	Daughters.Set(3);
	Daughters[0] =makeStaticData()->GetParticleID(d1);
	Daughters[1] =makeStaticData()->GetParticleID(d2);
	Daughters[2] =makeStaticData()->GetParticleID(d3);
    }
    else {
	if (!Next) Next = new PDecayChannel();
	Next->AddChannel(w,d1,d2,d3);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel:: AddChannel(Double_t w, char* d1, char* d2, char* d3, char* d4) {
    // The channel describes a decay into 4 daughter particles represented by
    // their names d1, d2, d3, d4 and branching ratio w.
    if (!Daughters.GetArray()) {
	Weight = w;
	NumOfDaughters = 4;
	Daughters.Set(4);;
	Daughters[0] =makeStaticData()->GetParticleID(d1);
	Daughters[1] =makeStaticData()->GetParticleID(d2);
	Daughters[2] =makeStaticData()->GetParticleID(d3);
	Daughters[3] =makeStaticData()->GetParticleID(d4);
    }
    else {
	if (!Next) Next = new PDecayChannel();
	Next->AddChannel(w,d1,d2,d3,d4);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void PDecayChannel:: AddChannel(Double_t w, Int_t nd, char** ld) {  
    // The channel describes a decay into nd daughter particles represented by
    // an array of their names ld and branching ratio w.
    if (!Daughters.GetArray()) {
	Weight = w;
	NumOfDaughters = nd;
	Daughters.Set(nd);;	
	for (Int_t i=0; i<nd; i++) {
	    Daughters[i] = makeStaticData()->GetParticleID(ld[i]);
	}
    }
    else {
	if (!Next) Next = new PDecayChannel();
	Next->AddChannel(w,nd,ld);
    }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
Double_t PDecayChannel:: GetWeight() {
    // Returns the branching ratio.
    return Weight;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
Int_t* PDecayChannel:: GetDaughters(int& number) {
    // The contents of 'number' is replaced by the number of daughter particles
    // and the pointer to the daughter array is returned.
    number = NumOfDaughters;
    return Daughters.GetArray();
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
PDecayChannel* PDecayChannel:: GetNext() {
    // Returns the 'Next' element, i.e. the succeeding channel.
    return Next;
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
PDecayChannel* PDecayChannel:: GetLast() {
    // Returns the 'last' element in the list.
  
    if (Next)
	return Next->GetLast();
    else return this;;
}
// -


// --------------------------------------------------------------------------
void PDecayChannel:: Print(const Option_t* delme)const {
    //
    // Prints the channel information.
    // 

    Int_t dalitz = 0;
  
    if (!NumOfDaughters) {
	cout << "Empty Channel." << endl;
	return;
    }

//    cout.width(8);
    cout << 100.*Weight << " %: " ;
    for (Int_t i=0; i<NumOfDaughters; i++) {
	if ( (dalitz!=2) && (i>0) ) cout << ", ";
	if (!(Daughters[i])) {
	    cout << "dalitz(";
	    dalitz = 2;
	}
	else {
	    cout << makeStaticData()->GetParticleName(Daughters[i]);
	    if (!(--dalitz)) cout << ")"; 
	}
    }
//    cout << endl;
    if (Next) {
	cout << ", ";
	Next->Print();
    }
    return;
}
// --------------------------------------------------------------------------

