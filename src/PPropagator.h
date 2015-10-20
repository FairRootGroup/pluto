// Author: H. Schuldes
// Written: 21.04.09

// Off-shell vector meson propagator

#ifndef _PPROPAGATOR_H_
#define _PPROPAGATOR_H_

#include "PChannelModel.h"
#include "PDynamicData.h"
#include "TComplex.h"

class PPropagator : public PChannelModel
{
public:

    PPropagator(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;
    
    using PChannelModel::GetAmplitude;
    
    TComplex GetAmplitude(Double_t *mass, Int_t *didx=NULL);
        
    void SetPID(int i) {
	// Overwrite pid if needed
	pid = i;  
    };
    
private:
    
    int pid; // Local PID to be overwritten
    ClassDef(PPropagator,0)  // Off-shell vector meson propagator

};

#endif // _PPROPAGATOR_H_
