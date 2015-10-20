// Author: B. Ramstein, I. Froehlich
// Written: 12.10.2008
// Revised: 

#include "PChannelModel.h"


class PDeltaDalitzFF : public PChannelModel  {
  
 public:

    using PDistribution::GetWeight;
    PDeltaDalitzFF(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char* delme=NULL) const;
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);
    
    void SetQED(Int_t i){
	//0 QED
	//1 VMD
	//-1 VMD mass-dependent width
	useQED=i;
    };
    void SetCC(Double_t a, Double_t b, Double_t c) {
	g_m2=a*a;
	g_e2=b*b;
	g_c2=c*c;
	rescale = 1;
    };
    

 private:

    Double_t facteur_rho(Double_t x);
    Double_t gk2(Double_t t);

    Double_t beta_prime;// rho meson coupling strength
    Double_t beta;// size parameter for N-Delta transition (fitted using space-like data)
    
    Double_t a2; // phase for space-like to time-like continuation (fitted to elastic time-like
    // form factor data) 

    //Masses:
    Double_t mass_pi,fm_pis ,mass_rho,gamma_rho,mass_delta,mass_n,mu_p;

    Double_t theta;
    

    //Coupling constants for QED version
    Double_t g_m2, g_e2, g_c2;
    Int_t    useQED, rescale,rhopid;
 

    ClassDef(PDeltaDalitzFF,0)  //Form factor model to be attached to Delta-Dalitz
};




