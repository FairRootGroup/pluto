/////////////////////////////////////////////////////////////////////////////
//
// PDecayChannel
//
//  Author:   Volker Hejny
//  Written:  27.08.99
//  Revised:  15.06.00 MK
//  Revised:  23.07.07 IF
//
// This class implements a list of decay modes (normally for one particle.
// Each instance is linked to a next one if there are more choices to 
// decay.
// -----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

#ifndef _PDECAYCHANNEL_H_
#define _PDECAYCHANNEL_H_

#include "TObject.h"
#include "TArrayI.h"
#include "TF1.h"
#include "TF2.h"
class PParticle;

using namespace std;
#include <iostream>



class PDecayChannel : public TObject {

 private:
    Double_t		Weight;			// branching ratio
    Int_t			NumOfDaughters;         // number of decay products
    TArrayI               Daughters;              // particle ids of decay products
    PDecayChannel*	Next;		        // pointer to an alternative channel
  
 public:
  
    PDecayChannel();
    ~PDecayChannel();
  
    // Add one decay channel by referencing an existing one.
    void AddChannel(PDecayChannel * n);

    // Add one decay channel by giving the weight and the particle ids. In
    // case of more than four particles one could also use a array of ids
    // together with its length.
    void AddChannel(Double_t w, Int_t d1); 
    void AddChannel(Double_t w, Int_t d1, Int_t d2);
    void AddChannel(Double_t w, Int_t d1, Int_t d2, Int_t d3);
    void AddChannel(Double_t w, Int_t d1, Int_t d2, Int_t d3, Int_t d4);
    void AddChannel(Double_t w, Int_t nd, Int_t* ld);
  
    // The same as above using PParticle* instead of ids.
    void AddChannel(Double_t w, PParticle* d1);
    void AddChannel(Double_t w, PParticle* d1, PParticle* d2);
    void AddChannel(Double_t w, PParticle* d1, PParticle* d2, PParticle* d3);
    void AddChannel(Double_t w, PParticle* d1, PParticle* d2, PParticle* d3, PParticle* d4);
    void AddChannel(Double_t w, Int_t nd, PParticle** ld);
  
    // The same as above using particle names instead of ids.
    void AddChannel(Double_t w, char* d1);
    void AddChannel(Double_t w, char* d1, char* d2);
    void AddChannel(Double_t w, char* d1, char* d2, char* d3);
    void AddChannel(Double_t w, char* d1, char* d2, char* d3, char* d4);
    void AddChannel(Double_t w, Int_t nd, char** ld);
  
    // Get the information out of the class ...
    Double_t GetWeight();
    Int_t* GetDaughters(int& n);
    PDecayChannel* GetNext();
    PDecayChannel* GetLast();

    void setAngleFunction(TF1 * f) {cout << "setAngleFunction is obsolete" << endl;};
    void setAngleFunction(TF2 * f) {cout << "setAngleFunction is obsolete" << endl;};
  
    // Unlink means to return the next decay channel and to set Next to NULL
    PDecayChannel* Unlink();
  
    // print decay channel information
    void Print(const Option_t* delme=NULL) const;
  
    ClassDef(PDecayChannel,0)//Pluto Decay Channel Class
	}; 
// end of PDecayChannel
						
#endif // _PDECAYCHANNEL_H_
