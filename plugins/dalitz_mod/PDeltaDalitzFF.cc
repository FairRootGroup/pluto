/////////////////////////////////////////////////////////////////////
// 
// Dalitz decay form factor models
//
// QED: Version with factors from M.I. Krivoruchenko, nucl-th/0104045
//      (Ref21)
//      The coupling constants can be changed via
//         SetCC(G_m,G_e,G_c);
// VMD: Calculation by F. Iachello and Q. Wan
//      Int. J. Mod. Phys.  A 20 (2005) 1846;
//      Q. Wan, Ph.D. Thesis, Yale University, New Haven, Connecticut (2007);
//      (Ref22)
//      Provided by B. Ramstein
//         SetCC(G_m,G_e,G_c); -> re-normalize VMD model such that
//                                it matches to the QED at pole mass
// 
//                                  Author:  I. Froehlich / B. Ramstein
/////////////////////////////////////////////////////////////////////



#include "PDeltaDalitzFF.h"
#include "PDynamicData.h"

PDistribution* PDeltaDalitzFF::Clone(const char*) const {
    //clone the object
    return new PDeltaDalitzFF((const PDeltaDalitzFF &)* this);
};

PDeltaDalitzFF::PDeltaDalitzFF(const Char_t *id ,const Char_t *de, Int_t key) : 
    PChannelModel(id, de,key) {
    //Constructor

    beta_prime = 0.004;
    beta       = 1.2748;
    a2         = 0.29;
    theta      = 53.;

    mass_pi    = makeStaticData()->GetParticleMass("pi+");
    fm_pis     = 4*mass_pi*mass_pi;
    //  mass_rho   = makeStaticData()->GetParticleMass("rho0");
    //  gamma_rho  = makeStaticData()->GetParticleTotalWidth("rho0");
      mass_rho   = 0.765;
      gamma_rho  = 0.112;
    mass_delta = makeStaticData()->GetParticleMass("D+");
    mass_n     = makeStaticData()->GetParticleMass("n");

    rhopid = makeStaticData()->GetParticleID("rho0");

    mu_p = 2.793;
    g_m2 = 3.2*3.2;
    g_e2 = 0.04*0.04;
    g_c2 = 0.2 *0.2;

    useQED = -2;
    rescale = 0;
} ;

Double_t PDeltaDalitzFF::gk2(Double_t t) {
    // intrinsic form_factor 
    Double_t denom2=1 + a2*a2*t*t + 2.*a2*t*cos(theta*TMath::Pi()/180.);
//cout <<"t= "<<t<< " g= "<<1-1./denom2<<endl;
    return 1/denom2;
}

Double_t PDeltaDalitzFF::facteur_rho(Double_t x) {
    // rho propagator term 

    //Iachello "old" version			
    if ((x < fm_pis) && (useQED != -2)) return mass_rho*mass_rho
					    /(mass_rho*mass_rho - x);
    //scaling to avoid "step"
    //  Double_t scale=0.91;
    Double_t scale=1.;
    Double_t alpha=1.;
    Double_t beta =0.;

    if (useQED == -1)
	gamma_rho  = makeDynamicData()->GetParticleTotalWidth(sqrt(x),rhopid);
    
    if  ((useQED == -2) &&  (x < fm_pis)) { //TODO: Make a more meaningfull flag name
	alpha = sqrt((fm_pis - x)/x)*(1.-2/TMath::Pi() *atan(sqrt((fm_pis - x)/x)));
	beta = 0.;
    } else {
	Double_t qnorm=x/fm_pis;
	Double_t logpart = (sqrt(x - fm_pis ) + sqrt(x))/(2*mass_pi);
	alpha= 2/TMath::Pi()*sqrt((x - fm_pis )/x)*log(logpart);
	beta = sqrt(pow(qnorm-1.,3)/qnorm);
    }
    
    //     Double_t logpart = (sqrt(x - fm_pis ) + sqrt(x))/(2*mass_pi);
    //     Double_t alpha= 2/TMath::Pi()*sqrt((x - fm_pis )/x)*log(logpart);
    Double_t denom_reel = mass_rho*mass_rho - x + (fm_pis -x)*gamma_rho*alpha/mass_pi;
    //    Double_t qnorm=x/fm_pis;
    //    Double_t beta = sqrt(pow(qnorm-1.,3)/qnorm);
    double denom_im = gamma_rho*4*mass_pi*beta;
    return scale*(mass_rho*mass_rho + 
	    8*gamma_rho*mass_pi/TMath::Pi())
	/sqrt(denom_reel*denom_reel+denom_im*denom_im);
}

Double_t PDeltaDalitzFF::GetWeight(Double_t *mass, Int_t *) {

    Double_t x=mass[0]*mass[0];
    Double_t m_d2=mass[1]*mass[1]; //mass of Delta sq.

    //q-dependent QED-Formfactor
    if (useQED == 1)
	return (g_m2 + 3 * g_e2 + (g_c2*x/(2*m_d2)));
    

    //VMD-Formfactor
    Double_t ff = 2.4637*gk2(-x)*(beta_prime + beta*facteur_rho(x));
    if (rescale) ff *= sqrt((4.76826/4.68927)*(g_m2/(3.2*3.2)));
    return ff*ff;

}


ClassImp(PDeltaDalitzFF)


