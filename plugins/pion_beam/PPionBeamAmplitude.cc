////////////////////////////////////////////////////////
//  PPionBeamAmplitude
//
//  This class contains Eqns. (15-18, 31-33) from 
//  Ref. 23
//
//                    Author:  H. Schuldes / IF
//                    Written: 21.04.09
//
////////////////////////////////////////////////////////



#include "PPionBeamAmplitude.h"



PDistribution* PPionBeamAmplitude::Clone(const char*) const
{
    //clone the object
    return new PPionBeamAmplitude((const PPionBeamAmplitude &)* this);
};

PPionBeamAmplitude::PPionBeamAmplitude(const Char_t *id, const Char_t *de, Int_t key) : 
    PChannelModel(id, de,key) {
    //Constructor
    dilepton=NULL;
    parent=NULL;
    RhoPropagator=NULL;
    OmPropagator=NULL;
    p_in=NULL;
    pion=NULL; 
    p_out=NULL;
    term=0;
    mode=0;
    monte_carlo=0;
};

Bool_t PPionBeamAmplitude::Init(void)
{   
    //looking for the mandatory particles 
    dilepton = GetParticle("dilepton");
    if (!dilepton) {
	Warning("Init","Dilepton not found");
	return kFALSE;
    }
  
    parent = GetParticle("parent");
    if (!parent) {
	Warning("Init","Parent not found");
	return kFALSE;
    }
 
    p_in= GetParticle("p_in");
    if (!p_in) {
	Warning("Init","p_in not found");
	return kFALSE;
    }

    pion= GetParticle("pion");
    if (!pion) {
	Warning("Init","pion not found");
	return kFALSE;
    }

    p_out= GetParticle("p_out");
    if (!p_out) {
	Warning("Init","p_out not found");
	return kFALSE;
    }
    
    if (pion->Is("pi+")) mode=1;

    return kTRUE;
}


Double_t PPionBeamAmplitude::GetWeight(void) {

    PParticle p_in_tmp(p_in);
    PParticle pion_tmp(pion);
    p_in_tmp.Boost(-parent->BoostVector());
    pion_tmp.Boost(-parent->BoostVector());


    p_bar_p = p_out->P();
    q_bar_p = dilepton->P();           //Monte-Carlo
    p_p = p_in_tmp.P();
    p_0 = p_in_tmp.E();  
    p_bar_0 = p_out->E();
    s = parent->M2();
  
    s= parent->M2();
    monte_carlo=1;
  
    Double_t dilepton_mass = dilepton->M();
    Double_t wert = GetWeight(&dilepton_mass);

    return wert;
}


Double_t   PPionBeamAmplitude::GetWeight(Double_t *mass, Int_t *didx) {

    Double_t q_bar_m = mass[0];
    Double_t q_bar_m2 = q_bar_m * q_bar_m;

    Double_t f_Rho = 0.036 ;
    Double_t f_Om = 0.011 ; //Coupling constants of the rho- and w-meson (real und positiv)


  

    Double_t M_Rho = makeStaticData()->GetParticleMass(41);
    Double_t M_Om = makeStaticData()->GetParticleMass(52);
    Double_t M_p = makeStaticData()->GetParticleMass(14);
    Double_t M_n = makeStaticData()->GetParticleMass(13);
    Double_t m_e = makeStaticData()->GetParticleMass(3);
    Double_t m_pi_minus = makeStaticData()->GetParticleMass(9);
    //    Double_t m_pi_plus = makeStaticData()->GetParticleMass(8);

    //  if (q_bar_m < 2* m_e) return 0.;
    if (q_bar_m < m_pi_minus) return 0.;
  
    if (q_bar_m > (sqrt(s) - M_p)) return 0.;

    if (!monte_carlo) {  //if Pluto is not in Monte-Carlo mode use kinematic eqn from Ref. 23
	p_p = sqrt(s) / 2. * sqrt(1. - 2. * (M_p*M_p  + m_pi_minus*m_pi_minus) / s + 
				  (M_p*M_p - m_pi_minus*m_pi_minus)*(M_p*M_p - m_pi_minus*m_pi_minus) / (s*s) );
	
	p_bar_p = sqrt(s) / 2. * sqrt(1. - 2. * (M_p*M_p  + q_bar_m2) / s + 
				      (M_p*M_p - q_bar_m2)*(M_p*M_p - q_bar_m2) / (s*s) );
	
		
	p_0 = sqrt(p_p*p_p + M_p*M_p);
	
	p_bar_0 = sqrt(p_bar_p*p_bar_p + M_n*M_n);
	
    }

    Double_t C_1_2_1_2_bar = ((p_0 + M_p) * (p_bar_0 + M_n)) / (4. * M_n * M_p) 
	* (1. + p_bar_p * p_bar_p / (3. * q_bar_m2));

    Double_t C_3_2_3_2_bar = (p_bar_0 + M_n) / (p_0 + M_p)  
	* ((p_p*p_p*p_p*p_p) / (2. * M_n * M_p) )
	* (1. + p_bar_p * p_bar_p / (3. * q_bar_m2));

    Double_t sqrt_s = sqrt(s);

    TComplex Rho_I_1_2_J_1_2(Graph_Rho_Re_T11->Eval(sqrt_s), Graph_Rho_Im_T11->Eval(sqrt_s));
    TComplex Rho_I_1_2_J_3_2(Graph_Rho_Re_T13->Eval(sqrt_s), Graph_Rho_Im_T13->Eval(sqrt_s));
    TComplex Rho_I_3_2_J_1_2(Graph_Rho_Re_T31->Eval(sqrt_s), Graph_Rho_Im_T31->Eval(sqrt_s));
    TComplex Rho_I_3_2_J_3_2(Graph_Rho_Re_T33->Eval(sqrt_s), Graph_Rho_Im_T33->Eval(sqrt_s));

    TComplex Rho_J_1_2,Rho_J_3_2;

    if (mode==0){
	Rho_J_1_2 = sqrt(2.) / 3. * (Rho_I_3_2_J_1_2 - Rho_I_1_2_J_1_2); 
	Rho_J_3_2 = sqrt(2.) / 3. * (Rho_I_3_2_J_3_2 - Rho_I_1_2_J_3_2);
	//Paper Eqn (15)
    } else {
	Rho_J_1_2 = sqrt(2.) / 3. * (Rho_I_1_2_J_1_2 - Rho_I_3_2_J_1_2); 
	Rho_J_3_2 = sqrt(2.) / 3. * (Rho_I_1_2_J_3_2 - Rho_I_3_2_J_3_2);
    }


    TComplex Om_J_1_2(Graph_Om_Re_T11->Eval(sqrt_s), Graph_Om_Im_T11->Eval(sqrt_s));
    TComplex Om_J_3_2(Graph_Om_Re_T13->Eval(sqrt_s), Graph_Om_Im_T13->Eval(sqrt_s));

    Om_J_1_2 = sqrt(2. / 3.) * Om_J_1_2;
    Om_J_3_2 = sqrt(2. / 3.) * Om_J_3_2;
    //Paper Eqn (16)

    TComplex Rho_J_1_2_Star = TComplex::Conjugate(Rho_J_1_2);
    TComplex Rho_J_3_2_Star = TComplex::Conjugate(Rho_J_3_2);
    TComplex Om_J_1_2_Star = TComplex::Conjugate(Om_J_1_2);
    TComplex Om_J_3_2_Star = TComplex::Conjugate(Om_J_3_2);


    if (RhoPropagator == NULL) {
	RhoPropagator = makeDynamicData()->GetParticleSecondaryModel("rho0","propagator");
	if (!RhoPropagator) {
	    Error("GetWeight","Rho propagator not found");
	    return 0.;
	}
    }
  
    if (OmPropagator == NULL) {
	OmPropagator = makeDynamicData()->GetParticleSecondaryModel("w","propagator");
	if (!OmPropagator) {
	    Error("GetWeight","Omega propagator not found");
	    return 0.;
	}
    }

    TComplex S_rho = RhoPropagator->GetAmplitude(&q_bar_m);
    TComplex S_w =   OmPropagator->GetAmplitude(&q_bar_m);
 
    TComplex S_rho_Star = TComplex::Conjugate(S_rho);
    TComplex S_w_Star = TComplex::Conjugate(S_w);

    TComplex Part_1 = (f_Rho*f_Rho / (M_Rho*M_Rho*M_Rho*M_Rho)) 
	* S_rho_Star * S_rho
	* (C_1_2_1_2_bar * Rho_J_1_2_Star * Rho_J_1_2
	   + C_3_2_3_2_bar * Rho_J_3_2_Star * Rho_J_3_2);

    TComplex Part_2 = (f_Rho*f_Om / (M_Rho*M_Rho*M_Om*M_Om)
		       * S_rho_Star * S_w
		       * (C_1_2_1_2_bar * Rho_J_1_2_Star * Om_J_1_2
			  + C_3_2_3_2_bar * Rho_J_3_2_Star * Om_J_3_2));

    TComplex Part_3 = (f_Om*f_Rho / (M_Rho*M_Rho*M_Om*M_Om)
		       * S_w_Star * S_rho
		       * (C_1_2_1_2_bar * Om_J_1_2_Star * Rho_J_1_2
			  + C_3_2_3_2_bar * Om_J_3_2_Star * Rho_J_3_2));

    TComplex Part_4 = (f_Om*f_Om / (M_Om*M_Om*M_Om*M_Om)) 
	* S_w_Star * S_w
	* (C_1_2_1_2_bar * Om_J_1_2_Star * Om_J_1_2
	   + C_3_2_3_2_bar * Om_J_3_2_Star * Om_J_3_2);


    Double_t Vorfaktor = (1./137.) /(6.* TMath::Pi()*TMath::Pi()) 
	* M_p*M_n / s 
	* (p_bar_p / p_p) * m_e*m_e * 
	(1. + q_bar_m2 / (2. * m_e*m_e) )
	* sqrt(1. - 4. * m_e*m_e / q_bar_m2);

 
    Double_t d_Sigma_d_q_bar2=0.;

    if (term == 0)
	d_Sigma_d_q_bar2 = Vorfaktor * (Part_1 + Part_2 + Part_3 + Part_4).Re(); 
    else if (term == 1)
	d_Sigma_d_q_bar2 = Vorfaktor * (Part_1).Re(); 
    else if (term == 2)
	d_Sigma_d_q_bar2 = Vorfaktor * (Part_4).Re(); 

    return 2. * q_bar_m * d_Sigma_d_q_bar2 / 2100. ;

};



ClassImp(PPionBeamAmplitude)
