/////////////////////////////////////////////////////////////////////
//
// Decay Model of Hadron -> Hadron(stable) + Hadron(stable)
// Such a case is e.g. the Delta->N+pi
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PHadronDecay.h"


ClassImp(PHadronDecay)

PHadronDecay::PHadronDecay()  {
} ;

PHadronDecay::PHadronDecay(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de,key) {

    if (is_channel<0)
	Warning("PHadronDecay","The model (%s) should be bound to CHANNELS only",de);
  
    //Get particles
    Int_t tid[11];
    tid[0]=10; 
    makeStaticData()->GetDecayModeByKey(key,tid); // retrieve current mode info

    //Parent ALWAYS important (also for the inherited classes)
    parent_id   = makeStaticData()->GetDecayParentByKey(key);
    parent_g0   = makeStaticData()->GetParticleTotalWidth(parent_id);
    parent_mass = makeStaticData()->GetParticleMass(parent_id);

    if (tid[0]!=2) 
	Warning("PHadronDecay","(%s):  Only 2 body decay",de);

    mass1=makeStaticData()->GetParticleMass(tid[1]);
    mass2=makeStaticData()->GetParticleMass(tid[2]);
    id1=tid[1];
    id2=tid[2];

    //Auto-set kinematic factors (for HadronWidth)
    //These values might be overwritten in inherited classes
    if (makeStaticData()->IsParticleMeson(parent_id)||
	PData::IsDelta(parent_id)) {
	use_fixed_delta=1;
	fixed_delta=0.09;
    } else use_fixed_delta=0;

    angular_l = (makeStaticData()->IsParticleMeson(parent_id)) ? 
	makeStaticData()->GetParticleSpin(parent_id)/2 : 
	PData::LPW(parent_id,id1,id2);
    if (!makeStaticData()->IsParticleMeson(parent_id)) cutoff_l = 1;
    else cutoff_l = 0; //cutoff=pow(cutoff,l+1); in HadronWidth

    if (makeStaticData()->IsParticleMeson(parent_id)) use_m0_over_m=1;
    else use_m0_over_m=0;

    if ((makeStaticData()->GetParticleID("D++") == parent_id) ||
	(makeStaticData()->GetParticleID("D+") == parent_id) ||
	(makeStaticData()->GetParticleID("D0") == parent_id) ||
	(makeStaticData()->GetParticleID("D-") == parent_id))
	use_m0_over_m=1;
    
    cutoff_version=0;

    version_flag |= VERSION_MASS_SAMPLING;  //Only one mass sampling in the PChannel

    w0=makeStaticData()->GetDecayBR(is_channel); //Weight used for mormalization
} ;

PDistribution* PHadronDecay::Clone(const char*delme) const {
    return new PHadronDecay((const PHadronDecay &)* this);
};

Bool_t PHadronDecay::Init(void) {
    //Init function called once for each PChannel
    

    return kTRUE;
}

Double_t PHadronDecay::EvalPar(const Double_t *x, const Double_t *params) {
    return Eval(x[0]);
}
 
Double_t PHadronDecay::Eval(Double_t x, Double_t y , Double_t z , Double_t t ) const
{
    Double_t res;
    Double_t mass[3];
    mass[0]=x; 
    mass[1]=mass1;
    mass[2]=mass2;

    if (draw_option==0) {
	return ((PChannelModel*)this)->GetWeight(mass);
	//return res;
    }
    if (draw_option==1) {
	((PChannelModel*)this)->GetWidth(x,&res);
	return res;
    }
    if (draw_option==2) {
	((PChannelModel*)this)->GetBR(x,&res);
	return res;
    }
    return 0;
}

int PHadronDecay::GetDepth(int i) {
    
    makeStaticData()->SetDecayEmin(is_channel, mass1+mass2);
    return 0; //2 stable products -> depth is 0
}

Bool_t PHadronDecay::SampleMass(void) {
    //Mass-sampling wrapper

    return kTRUE;
};

Bool_t PHadronDecay::SampleMass(Double_t *mass, Int_t *didx) {
    //Not much to do here...
    //Since we have 2 stable products, for completeness
    //we reset the masses to the nominal value
    mass[1]=mass1;
    mass[2]=mass2;
    return kTRUE;
};


// Bool_t PHadronDecay::GetBR(Double_t mass, Double_t *br, Double_t totalwidth) {
//     //Calculates the mass-dependent br
//     //the totalwidth may be set by the user, otherwise we use the static one
//     //We take the 2body-ps into account

   
//     if (!GetWidth(mass,&width)) return kFALSE;
//     *br = width/totalwidth;

//     double sc_pole=

//     return kTRUE; 
	
// }


Bool_t PHadronDecay::GetWidth(Double_t mass, Double_t *width, Int_t didx) {

    if (makeStaticData()->GetPWidx(is_channel)==-1) {
	*width = parent_g0;
	return kFALSE;
    }
    *width = w0*HadronWidth(mass, mass1, mass2);
    return kTRUE;

}



// Comment (IF): This function was used in Width1, but it seems to be useless
// in the case of 2 stable products
// It just makes checks which are repeated in HWidth (now: HadronWidth)

// double PHadronDecay::HW(const double & ecm, const int & id, const int & ia=0, const int & ib=0) {
//     // for stable products

//     const double mproton=makeStaticData()->GetParticleMass("p"), 
// 	m2pi0=2.*makeStaticData()->GetParticleMass("pi0");
//     double m0=makeStaticData()->GetParticleMass(id), ma, mb, ms;
    
//     if (!(ia+ib)) {             // flag identifying N + 2-pion s-wave production
//       ma=mproton;               // in this case treat as two-body decay: 
//       mb=m2pi0;                 // 1st product is N, 2nd quasi-product is di-pion
//     } else {                    // two-hadron decay
//       ma=makeStaticData()->GetParticleMass(ia);  // 1st hadron mass
//       mb=makeStaticData()->GetParticleMass(ib);  // 2nd hadron mass
//     }
//     ms=ma+mb;                   // sum of product masses

//     return (ecm<ms||m0<ms) ? 0. : HWidth(id,ia,ib,ecm,ma,mb);
//   }


double PHadronDecay:: HadronWidth(const double & m, const double & ma, const double & mb) {
    // Kinematic width (+ phase space + cutoff) for decays of type hadron -> hadron + hadron
    // Arguments: masses (GeV/c**2)
    // See Ref 10 Eqs. (15-17).
    //
    // Uses data members:
    // * parent_mass
    // * parent_g0
    
    double m0=parent_mass, 
	ms=ma+mb, 
	md=m0-ms;
//    if (m<=ms||md<=0.) return 0.;        // kinematically inaccessible
    if (m<=ms) return 0.;


    double qr2=PKinematics::pcms2(m0,ma,mb),
	q2=PKinematics::pcms2(m,ma,mb);
    if (qr2==0) //Below threshold production
	qr2=1;
    
    double q=sqrt(q2/qr2);

//    return q*parent_g0;

    if (!id1) return parent_g0*q*m0/m;  // ia=ib=0 is two pions coupled to s-wave

    double delta2=(use_fixed_delta) ? fixed_delta : md*md+0.25*parent_g0*parent_g0;
    double cutoff = (qr2 + delta2)/(q2 + delta2);

    if (angular_l) {
	q=pow(q,2*angular_l+1);
	if (cutoff_l) cutoff=pow(cutoff,angular_l+1);
    }

    if (cutoff_version==1) {
	if (angular_l)
	    cutoff=(1.2)/(1+0.2*pow(q,2*angular_l));
	else
	    cutoff=(1.2)/(1+0);
	//cout << angular_l << endl;
    } else if (cutoff_version==2) {
	cutoff=1;
	cout << "no cutoff " << is_channel << endl;
	cout << "mass " << m << " res: " << q << endl;
    }

    q*=parent_g0*cutoff;

//    cout << q << endl;
    return (use_m0_over_m) ? q*pow((m0/m),use_m0_over_m) : q;
}

ClassImp(PHadronDecay)
