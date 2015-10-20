// Author: H. Schuldes
// Written: 29.05.09


#ifndef _PPIONBEAMAMPLITUDE_H_
#define _PPIONBEAMAMPLITUDE_H_

#include "PChannelModel.h"
#include "PDynamicData.h"
#include "TComplex.h"
#include "TGraph.h"

class PPionBeamAmplitude : public PChannelModel
{
 public:

    PPionBeamAmplitude(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;
    
    using PChannelModel::GetWeight;
  
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);
    Bool_t   Init();
    Double_t GetWeight(void);

    void SetGraph_Rho_Re_T11(TGraph * f) {
	Graph_Rho_Re_T11=f;
    };
 
    void SetGraph_Rho_Im_T11(TGraph * f) {
	Graph_Rho_Im_T11=f;
    };   

    void SetGraph_Rho_Re_T13(TGraph * f) {
	Graph_Rho_Re_T13=f;
    };

    void SetGraph_Rho_Im_T13(TGraph * f) {
	Graph_Rho_Im_T13=f;
    };

    void SetGraph_Rho_Re_T31(TGraph * f) {
	Graph_Rho_Re_T31=f;
    };

    void SetGraph_Rho_Im_T31(TGraph * f) {
	Graph_Rho_Im_T31=f;
    };

    void SetGraph_Rho_Re_T33(TGraph * f) {
	Graph_Rho_Re_T33=f;
    };

    void SetGraph_Rho_Im_T33(TGraph * f) {
	Graph_Rho_Im_T33=f;
    };

    void SetGraph_Om_Re_T11(TGraph * f) {
	Graph_Om_Re_T11=f;
    };

    void SetGraph_Om_Im_T11(TGraph * f) {
	Graph_Om_Im_T11=f;
    };

    void SetGraph_Om_Re_T13(TGraph * f) {
	Graph_Om_Re_T13=f;
    };

    void SetGraph_Om_Im_T13(TGraph * f) {
	Graph_Om_Im_T13=f;
    };

    void SetS(double f) {
	s=f;
    };

    void SetTerm(int f) {
	//terms: 0=all
	//       1=rho
	//       2=omega
	term=f;
    };

    void SetMode(int f) {
	//mode: 0=pi_minus
	//      1=pi_plus  
	mode=f;
    };


 private:
    
    TGraph * Graph_Rho_Re_T11;
    TGraph * Graph_Rho_Im_T11;
    TGraph * Graph_Rho_Re_T13;
    TGraph * Graph_Rho_Im_T13;
    TGraph * Graph_Rho_Re_T31;
    TGraph * Graph_Rho_Im_T31;
    TGraph * Graph_Rho_Re_T33;
    TGraph * Graph_Rho_Im_T33;
    TGraph * Graph_Om_Re_T11;
    TGraph * Graph_Om_Im_T11;
    TGraph * Graph_Om_Re_T13;
    TGraph * Graph_Om_Im_T13;
    
    PChannelModel *RhoPropagator;
    PChannelModel *OmPropagator;
    
    PParticle *dilepton, *parent, *p_in, *pion, *p_out;

    Double_t p_bar_p; //Momentum of outgoing Neutron
    Double_t q_bar_p; //Momentum of Dilepton
    
    Double_t p_0;     //Energy of incoming Proton    
    Double_t p_bar_0; //Energy of Neutron
    Double_t p_p;     //Momentum of Proton
    Double_t s;       //Mass of parent squared
    Int_t term, monte_carlo;
    Int_t mode;

    Double_t f_Rho, f_Om; //Coupling constants of the rho- and w-meson (real und positiv)
    

    Double_t M_Rho, M_Om, M_p, M_n, m_e, m_pi_minus, m_pi_plus; //some constants
     
    ClassDef(PPionBeamAmplitude,0)  // Pion beam amplitudes

};

#endif // _PPIONBEAMAMPLITUDE_H_
