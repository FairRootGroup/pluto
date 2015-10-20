/////////////////////////////////////////////////////////////////////
//
//  N + N --> N + Delta, Ref 1
//
//                                  Author: Kagarlis  
//                                  Reimplemented I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "PDeltaAngularDistribution.h"

ClassImp(PDeltaAngularDistribution)

PDeltaAngularDistribution::PDeltaAngularDistribution() {

} ;

PDeltaAngularDistribution::PDeltaAngularDistribution(const Char_t *id, const Char_t *de) :
    PAngularDistribution(id, de) {

    beam = NULL;
    target = NULL;
    use_term = 1+2+4; //By default all terms
} ;

PDistribution* PDeltaAngularDistribution::Clone(const char*delme) const {
    return new PDeltaAngularDistribution((const PDeltaAngularDistribution &)* this);
};


Bool_t PDeltaAngularDistribution::IsValid(void) {
    
    if (!direct_sampling_done) {
	Fatal("IsValid","Sampling not finished");
    }

    return kTRUE;


};

Bool_t PDeltaAngularDistribution::Init(void) {
    if (!PAngularDistribution:: Init()) return kFALSE;

    //In addition we need the beam and target
    //done already in PAngularDistribution

    if (!beam || !target) {
	Warning("Init","beam or target not found");
	return kFALSE;
    }

    return kTRUE;
}

Bool_t PDeltaAngularDistribution::Prepare(void) {
    getNN_DeltaN_param();
    return kTRUE;
}

void PDeltaAngularDistribution::getNN_DeltaN_param() {
    // parameters of the Delta cm production angle in N+N->N+Delta
    // Ref. 1

    const double L1=0.3969, L2=0.36;     // NN & NDelta coupling constant^2
    const double fs_pi=2.202, 
	mpi=makeStaticData()->GetParticleMass("pi0"), g_pi=0.6, 
	fmg=g_pi*fs_pi/mpi, 
	mp=makeStaticData()->GetParticleMass("p"), mp2=mp*mp, mp4=mp2*mp2;
    double I2=(beam->Vect4())*(target->Vect4());// product of Lorentz vectors

//  if(pdNNN>0) I2 /= 2.;                 // for p+d and d+p, half I2
// Commented out, since treated explicetly


    I2=I2*I2-mp4;                         // kinematical factor for cross section
    if (I2<=0.) {
	anisotropy=0;
	return;
    }
    anisotropy=1;
    // factors for N + N --> N + Delta cross section
    lambda2 = ((beam->Vect()).Mag()<3.)?L1:L2;
    prefac = fmg*fmg/(4.*64.*TMath::Pi()*I2);

}



double PDeltaAngularDistribution::ds_dt(double cos_th_cm) {
    // ds/dt(cos_th_cm) in the cm for N+N->N+Delta. With the mass of Delta
    // sampled independently in PData, the simulated events are distributed
    // as ds/dOmega (normalized). (see Ref. 1)

    const double mpi=makeStaticData()->GetParticleMass("pi0"), 
	mpi2=mpi*mpi, 
	mp=makeStaticData()->GetParticleMass("p"),
	mp2=mp*mp, mp4=mp2*mp2, 
	mdelta=makeStaticData()->GetParticleMass("D0");

    double m=mres, md2=m*m, md4=md2*md2, mdmn=m-mp, mdn=m+mp,
	mdmnm=mdmn*mdn, mdn2=mdn*mdn, mdn4=mdn2*mdn2, mdmn2=mdmn*mdmn,
	pf=prefac/md2;

    // Mandelstam invariant t and u:
    double t=parent->InvariantT(mp,m,-cos_th_cm), // cos_th_cm sign inversed, Tingting Liu, 2010-04-13
	u=parent->InvariantT(m,mp,cos_th_cm),  // cos_th_cm sign inversed, Tingting Liu, 2010-04-13
	tu=t*u, tpu=t+u;

    // Form factors for t and u channels
    double f_t=(lambda2-mpi2)/(lambda2-t);       // form factor F(t) eq. (4)
    f_t*=f_t/(t-mpi2);
    double f_u=(lambda2-mpi2)/(lambda2-u);       // same for u channel F(u)
    f_u*=f_u/(u-mpi2);

    // Off-shell corrections of Eq. (17) of Ref 1,
    // but for Moniz Delta-width parametrization see Z. Phys. A356 (1997) 421
    double qt=(PKinematics::pcmt(mdelta,t)+.09)/(PKinematics::pcmt(m,t)+.09),
	qu=(PKinematics::pcmt(mdelta,u)+.09)/(PKinematics::pcmt(m,u)+.09);
    f_t *= qt*qt;                                // Eqs. (17)
    f_u *= qu*qu;
  
    // the matrix elements for N+N->N+Delta, Eqs. (7-8) in Ref 1
    double M_t2 = pf/3. * f_t * f_t * t * (t-mdmn2) * (t-mdn2) * (t-mdn2); // t-channel
    double M_u2 = pf/3. * f_u * f_u * u * (u-mdmn2) * (u-mdn2) * (u-mdn2); // u-channel
    double M_tu = pf/2. * f_t * f_u *
	( (tu + mdmnm * tpu - md4 + mp4) * (tu + mp * mdn * mdmnm)
	  - (tu - mdn2 * tpu + mdn4) * (tu - mp * mdmn * mdmnm)/3. ); // exchange term
    // fixed according to S. Teis Ph.D.  (Giessen 1996) 
    double M2 = 0;
    if (use_term & 0x1) {M2+=M_t2;primary->SetValue(T_MATRIX ,M_t2);}
    if (use_term & 0x2) {M2+=M_tu;primary->SetValue(TU_MATRIX ,M_tu);}
    if (use_term & 0x4) {M2+=M_u2;primary->SetValue(U_MATRIX ,M_u2);}
  
    return M2;
}

double PDeltaAngularDistribution::SamplePolarAngle(double r) {

    double f=-1, tf=-1, r2, c0=-1, r1=r, a=-1, b=-1, area=-1, x1, delta;
    int fl=0 ;

    if (!anisotropy) return 2*r-1;

    mres=primary->M();

    // test-function parameters for resonances
    a = 1.01*ds_dt(0.);                 // constant
    b = 1.2*(1.01*ds_dt(.995) - a);     // slope  (ds_dt(0) droops, use 0.995)
    area = b*(2. + b/a)/a;              // area x b/a^2
  
 again:
  
    fl=1;
    x1=(1-2*(r1>.5))*a/b;
    delta=sqrt(1+area*fabs(1-2*r1));
    c0=x1*(1-delta);
    f  = ds_dt(c0);
    tf = a+b*fabs(c0);
//  printf("********a=%f b=%f: c0=%f tf=%f f=%f\n",a,b,c0,tf,f);
    r2=PUtils::sampleFlat();
    if (r2>f/tf) {                        // reject
	r1=PUtils::sampleFlat();
	goto again;
    }
    if (tf<f) {                           // diagnostics
//	Warning("SamplePolarAngle","algorithm failed %d",fl);
//	Warning("SamplePolarAngle","a=%f b=%f: c0=%f tf=%f f=%f mres=%f mpar=%f",a,b,c0,tf,f,mres,parent->M());
	r1=PUtils::sampleFlat();
	goto again;
    }
    return c0;                            // accept
  

}

void PDeltaAngularDistribution::Print(const Option_t* delme) const {

    PAngularDistribution::Print();
}


