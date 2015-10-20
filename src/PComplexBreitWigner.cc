/////////////////////////////////////////////////////////////////////
//
// Resonance mass sampling from a 
// complex breit-wigner propagator
// for each daughter decay channel, a list of
// interference channel can be set up
// This includes the complex amplitude for the given
// decay channel as well as the interference terms
// The interference terms can be other particle models,
// decay models and even exchange terms itself
//
// Such, by removing the standard propagator (setting the amplitude to 0)
// and adding other exchange terms, the production
// of a single particle (with known target daugters) can be rewritten
// in forms of koherent feynman sums
//
//                                  Author: Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PComplexBreitWigner.h"


ClassImp(PComplexBreitWigner)

PComplexBreitWigner::PComplexBreitWigner()  {

} ;

PComplexBreitWigner::PComplexBreitWigner(const Char_t *id, const Char_t *de, Int_t key) :
    PBreitWigner(id, de, key) {

    readModesDone=0;
    readModelsDone=0;
    updateAmplitude=0;
    mr=makeStaticData()->GetParticleMassByKey(key);
};

PDistribution* PComplexBreitWigner::Clone(const char*delme) const {
    return new PComplexBreitWigner((const PComplexBreitWigner &)* this);
};

Bool_t PComplexBreitWigner::SampleMass(Double_t *mass, Int_t *didx) {
    // Resonance mass sampling from a relativistic
    // Breit-Wigner 

    if (didx) 
	didx_option=didx[0];
    else
	didx_option=-1;
    mass[0] = this->GetRandom();
    return kTRUE;
}

Double_t PComplexBreitWigner::GetWeight(Double_t *mass, Int_t *didx) {
    return GetAmplitude(mass, didx).Rho2();
}


void PComplexBreitWigner::ReadModes(void) {
    if (readModesDone) return;

   //on the first call, loop over decay modes
    Int_t listkey=-1;

    num_decaychannels=0;
    while (makeDataBase()->MakeListIterator(primary_key, "pnmodes", 
					    "link", &listkey)) {
	
	int_index[num_decaychannels]=makeStaticData()->GetDecayIdxByKey(listkey);
	num_terms[num_decaychannels]=0;
	//By default loacl amplitude is (1,0)
	//Do not forget term "0" is always the local term
	int_phase[num_decaychannels][0]=0;
	int_ampl[num_decaychannels][0]=1;
	
	num_decaychannels++;

	if (num_decaychannels==COMPLEX_MAX_DECAYCHANNELS) {
	    Warning("ReadModes","COMPLEX_MAX_DECAYCHANNELS reached");
	    listkey=-1;
	}
    } //end iterator
    readModesDone=1;
}

void PComplexBreitWigner::AddAmplitude(int idx,double ampl,double phase) {
    ReadModes();

    //loop over modes and try to find the idx

    int found_idx=-1;
    for (int i=0;i<num_decaychannels;i++) 
	if (int_index[i]==idx) found_idx=i;
    if (found_idx<0) {
	Warning("AddAmplitude:","idx %i is not a valid decay index for particle %i",idx,is_pid);
	return;
    }

    int_phase[found_idx][0]=phase;
    int_ampl[found_idx][0]=ampl;

}

void PComplexBreitWigner::AddInterference(int idx,int key,int didx,double ampl,double phase) {
    // Adds a complex interference to the decay index "idx"
    // The "idx" must be a valid decay for the particle which is assigned to the "PComplexBreitWigner"
    // model. "key" is the data base entry wihc is used for the mixing (e.g. another complex BW)
    // "didx" is the target decay index used with the "key" model to get the correct final state
    // "ampl" and "phase" is the complex amplitude given in polar coordinates
    ReadModes();

    //loop over modes and try to find the idx

    int found_idx=-1;
    for (int i=0;i<num_decaychannels;i++) 
	if (int_index[i]==idx) found_idx=i;
    if (found_idx<0) {
	Warning("AddInterference","idx %i is not a valid decay index for particle %i",idx,is_pid);
	return;
    }

    if (num_terms[found_idx]==COMPLEX_MAX_TERMS){
	Warning("AddInterference","COMPLEX_MAX_TERMS reached");
	return;
    }

    num_terms[found_idx]++;

    int_didx[found_idx][num_terms[found_idx]]=didx;
    int_key[found_idx][num_terms[found_idx]]=key;
    int_phase[found_idx][num_terms[found_idx]]=phase;
    int_ampl[found_idx][num_terms[found_idx]]=ampl;

};

Double_t PComplexBreitWigner::EvalPar(const Double_t *x, const Double_t *params) {
    if (params) {
	draw_option=(int)params[0];
	didx_option=(int)params[1];

	if (updateAmplitude) {
	    int found_idx=-1;	    
	    //now look for didx
	    for (int i=0;i<num_decaychannels;i++) 
		if (int_index[i]==didx_option) found_idx=i;
	    if (found_idx>=0) {int_phase[found_idx][1]=params[3];
	    int_ampl[found_idx][1]=params[2];}
	}
    }
    return Eval(x[0]);
}
 
Double_t PComplexBreitWigner::Eval(Double_t x, Double_t y , Double_t z , Double_t t ) const
{
    Double_t res;
    Double_t my_x[11];
    Int_t my_didx_option[11];
    my_x[0]=x;
    my_didx_option[0]=didx_option;

    if (draw_option==0) {
	return ((PChannelModel*)this)->GetWeight(my_x,my_didx_option);
	//return res;
    }
    if (draw_option==1) {
	((PChannelModel*)this)->GetWidth(x,&res,didx_option);
	return res;
    }
    if (draw_option==2) {
	((PChannelModel*)this)->GetBR(x,&res);
	return res;
    }
    if (draw_option==3) {
	return ((PChannelModel*)this)->GetAmplitude(my_x,my_didx_option).Re();
    }
    if (draw_option==4) {
	return ((PChannelModel*)this)->GetAmplitude(my_x,my_didx_option).Im();
    }
    return 0;
}

void PComplexBreitWigner::ReadModels(void) {
    if (readModelsDone) return;

    //loop over decay modes
    for (int i=0;i<num_decaychannels;i++) {
	//Now loop over possible interference terms
	for (int j=1;j<=num_terms[i];j++) {
	    TObject *t_result;
	    if (makeDataBase()->GetParamTObj (int_key[i][j] ,"model", &t_result))
		p[i][j]=(PChannelModel*)t_result;
	    else {
		p[i][j]=NULL;	
		Warning("ReadModels","No term model found: this will not work");
	    }
	    p[i][j]=(PChannelModel*)t_result;
	}
    }
    readModelsDone=1;
}


TComplex PComplexBreitWigner::GetAmplitudeLocal(Double_t *mass, Int_t num) {
    
    int didx_local=int_index[num];
    double m=mass[0];
    double m2=m*m, 	
	mm=mr*mr-m2;

    global_weight_scaling = makeDynamicData()->GetParticleScalingFactor(is_pid);

    double wmt = 1.;

    if (!GetWidth(m,&width)) return -1;
 

    double partial_width;
    if (didx_local>=0) {
	if (!GetWidth(m,&partial_width,didx_local)) return -1;
	if (partial_width > width*1.1) {
	}
    } else partial_width=width;

    TComplex denom(mm,m*width);
    
    TComplex numerator(m*sqrt(partial_width*global_weight_scaling*wmt),0);
    
    numerator/=denom;
    TComplex ampl(int_ampl[num][0],int_phase[num][0],kTRUE); //local amplitude
    numerator*=ampl;
    
    
    //Now loop over possible interference terms

    TComplex interference(0,0);
    for (int i=1;i<=num_terms[num];i++) {
	TComplex interference_term(int_ampl[num][i],int_phase[num][i],kTRUE);

	if (p[num][i]) {

	    interference_term*=p[num][i]->GetAmplitude(mass,&int_didx[num][i]);
	} else {
	    Warning("GetAmplitudeLocal","No model found");
	}
	interference+=interference_term;
    }

    if (num_terms[num])
	return numerator+interference;
    return numerator;

}


TComplex PComplexBreitWigner::GetAmplitude(Double_t *mass, Int_t *didx) {
// relativistic Breit-Wigner distribution function for particle "key"
    ReadModes();
    ReadModels();
    double weight=0;
    int didx_local=-1;
    if (didx) didx_local=didx[0];

    //loop over decay modes
    for (int i=0;i<num_decaychannels;i++) {
	if (didx_local<0)
	    weight+=GetAmplitudeLocal(mass,i).Rho();
	else if (didx_local==int_index[i])
	    return GetAmplitudeLocal(mass,i);
    }
    TComplex ampl(weight,0);
    return ampl;

 }


void PComplexBreitWigner::Print(const Option_t* delme) const {
    //loop over decay modes
    for (int i=0;i<num_decaychannels;i++) {
	//Now loop over possible interference terms
	cout << "Num:" << i << " Didx: " << int_index[i] << endl; 
	for (int j=1;j<=num_terms[i];j++) {
	    cout << " Term: " << j << endl; 
	}

    }
}
