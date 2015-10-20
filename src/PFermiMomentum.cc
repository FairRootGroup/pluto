/////////////////////////////////////////////////////////////////////
//
// Quasi-free scattering and construction of virtual quasi-beam-particle
//
//                                  Author: R. Holzmann
//                                  Reimplemented I. Froehlich
// Ref 1: Benz et al. NP B65 (1973) 158
//        (p + d Fermi kinematics, implemented by R. Holzmann)
//
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PFermiMomentum.h"
#include "PDynamicData.h"

ClassImp(PFermiMomentum)

    PFermiMomentum::PFermiMomentum() {

} ;

PFermiMomentum::PFermiMomentum(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de, key) {

    beam = NULL;
    target = NULL;
    spectator = NULL;
    parent = NULL;
    p1 = NULL;
    p2 = NULL;
    composite = NULL;
    SetNpx(500);
    fXmin      = 0.;
    fXmax      = 1.;
    mesh=NULL;

    spline=kFALSE;
    g_spline=NULL;
    mom=NULL;

    didx_composite=-1;
    composite_model=NULL;
    tcross_key=makeStaticData()->MakeDirectoryEntry("modeldef",NMODEL_NAME,LMODEL_NAME,"tcross");
    relative_warning=0;

} ;

PDistribution* PFermiMomentum::Clone(const char*delme) const {
    return new PFermiMomentum((const PFermiMomentum &)* this);
};


PFermiMomentum::~PFermiMomentum() {

    if (mesh) delete mesh;
}

double PFermiMomentum:: SampleFermi(double & px, double & py, double & pz) {
    // sample Deuteron Fermi momentum in GeV/c
  
    double p = this->GetRandom();
    while (p>0.3)
	p = this->GetRandom();

    // cout << p << endl;

    double theta = acos(1.-2.*PUtils::sampleFlat());
    double phi = 2.*TMath::Pi()*PUtils::sampleFlat();
    double sth=sin(theta);
    px = p*cos(phi)*sth;
    py = p*sin(phi)*sth;
    pz = p*cos(theta);

    return p;
}


Double_t PFermiMomentum:: GetWeight(Double_t *mass, Int_t * ) {
    //Use GetWeight for multi-dimensional sampling
    //mass[0]: Momentum of nucleon
    //mass[0]: cos Theta (polar angle)

    double p = mass[0];
    double theta = acos(mass[1]);
    double phi = 2.*TMath::Pi()*PUtils::sampleFlat();
    double sth=sin(theta);
    px = p*cos(phi)*sth;
    py = p*sin(phi)*sth;
    pz = p*cos(theta);

    Double_t massS, eS, eP, ptot, t=-1., mdeut;

    ptot = sqrt(px*px+py*py+pz*pz);

    massS = spectator->M();                        // mass of spectator nucleon
    eS = sqrt(ptot*ptot + massS*massS);            // spectator total energy in deuteron c.m.
    mdeut = makeStaticData()->GetParticleMass("d");
    t = pow(mdeut-massS,2) - 2.*mdeut*(eS-massS);  // off-shell mass**2 of participant
    eP = sqrt(ptot*ptot + t);         // participant total energy	
    TLorentzVector my_part(px,py,pz,eP);

    TLorentzVector my_q = my_part + my_beam;
    Double_t my_sqrt = my_q.M();



    Double_t w = Eval(p);
    if (composite_model && (!(composite_model->GetVersionFlag()) & VERSION_WEIGHTING)) {
	w *= composite_model->GetWeight(&my_sqrt);
    }

    // if (debug_print) cout << "my_sqrt: " << my_sqrt << "w: " << composite_model->GetWeight(&my_sqrt) << endl;

    return w;
    
}


double PFermiMomentum:: EvalPar(const double* x, const double* par) {
    return Eval(x[0]);
}

Double_t PFermiMomentum::Eval(Double_t x, Double_t y , Double_t z , Double_t t ) const {
    // Deuteron wave function in p-space (Ref 6)

    double p = x;
    double r = p/0.197;

    if (mom) {
	Double_t x = mom->Eval(r,g_spline,"");
	if (x>0) return r*r*x;
	return 0;
    }

    const double sqrtpi2 = 0.79788456;
    const double alpha = 0.23162461;

    double m[13], m2[13];
    for(int i=0;i<13;i++) {
	m[i] = alpha + i;
	m2[i] = m[i]*m[i];
    }
  
    double c[13] = {0.88688076e0, -0.34717093e0, -0.30502380e1, 0.56207766e2,
		    -0.74957334e3,  0.53365279e4, -0.22706863e5, 0.60434469e5,
		    -0.10292058e6,  0.11223357e6, -0.75925226e5, 0.29059715e5,
		    0.};
  
    double d[13] = {0.23135193e-1, -0.85604572e0, 0.56068193e1, -0.69462922e2,
		    0.41631118e3,  -0.12546621e4, 0.12387830e4,  0.33739172e4,
		    -0.13041151e5,   0.19512524e5, 0., 0., 0.}; 
  
    c[12] = 0.;
    for(int i=0;i<12;i++) c[12] -= c[i];   // normalize c[] properly
  
    int n = 12, n1 = 11, n2 = 10;
  
    double sum1 = 0.;
    double sum2 = 0.;
    double sum3 = 0.;
    double rtemp;
    for(int i=0;i<5;i++) {
	rtemp = d[i*2]/m2[i*2] + d[i*2+1]/m2[i*2+1];
	sum1 = sum1 + rtemp;
	rtemp = d[i*2] + d[i*2+1];
	sum2 = sum2 + rtemp;
	rtemp = d[i*2]*m2[i*2] + d[i*2+1]*m2[i*2+1];
	sum3 = sum3 + rtemp;
    }
  
    for(int i=0;i<3;i++) {                 // normalize d[] properly
	d[n2] = -m2[n1]*m2[n]*sum1 + (m2[n1]+m2[n])*sum2 - sum3;
	d[n2] = d[n2] * m2[n2]/(m2[n]-m2[n2])/(m2[n1]-m2[n2]);
	int cycle = n2;
	n2 = n1;
	n1 = n;
	n = cycle;
    }
  
    double U = 0.;
    double W = 0.;
    for(int i=0;i<13;i++) {
	U += c[i]/(r*r + m2[i]);
	W += d[i]/(r*r + m2[i]);
    }
    U = sqrtpi2 * U;    // s wave
    W = sqrtpi2 * W;    // d wave
  
    return r*r*(U*U + W*W);  // total probability to have momentum p
}


Bool_t PFermiMomentum::IsValid(void) {
    
   
    return kTRUE;


};

Bool_t PFermiMomentum:: Init(void) {

    num_of_sampledevents=num_of_realevents=0; //BUGBUG: Can cause trouble in BulkDecay

    beam = GetParticle("beam");
    target = GetParticle("target");
    parent=GetParticle("parent");

    if (!beam || !target) {
	Warning("Init","<%s> beam or target not found",GetIdentifier());
	return kFALSE;
    }

    spectator= GetParticle("spectator");
    p1=GetParticle("p1");
    p2=GetParticle("p2");
    //BUGBUG: Check presence
    composite= GetParticle("composite");
    if (!composite) {
	Warning("Init","<%s> composite not found",GetIdentifier());
	return kFALSE;
    }

    return kTRUE;
}

int PFermiMomentum::GetDepth(int i) {
    return 0;
}

void PFermiMomentum::SubPrint(Int_t opt) const {
    //Print sub-models    
    if (composite_model) {cout << " {";cout << composite_model->GetIdentifier()<<"}";}
}


Bool_t PFermiMomentum:: Prepare(void) {

    if (composite) {
	didx_composite = composite->GetDecayModeIndex(1);
	// cout << didx_composite << endl;
	if (didx_composite>0) {
	    composite_model=
		makeDynamicData()->GetDecayModelByKey(
		    makeStaticData()->GetDecayKey(didx_composite),tcross_key);
	}
	
	//   if (composite_model) composite_model->Print();
    } 

    return kTRUE;
}

Bool_t PFermiMomentum:: SampleMass(void) {

    
    Int_t beam_breakup=0;

    //This is the d+N reaction:
    if ((beam->ID() == makeStaticData()->GetParticleID("d")) && (
	    target->ID() != makeStaticData()->GetParticleID("d"))) {
	beam_breakup=1;
    } 
    else if ((beam->ID() == makeStaticData()->GetParticleID("d")) && (
	    target->ID() == makeStaticData()->GetParticleID("d"))) { 
	//d+d_coherent reaction
	Double_t is_breakup = 0;
	if (beam->GetValue(IS_BREAKUP, &is_breakup)) //value used
	    beam_breakup=1;
	else if (!target->GetValue(IS_BREAKUP, &is_breakup)) {
	    //randomize if value not used
	    beam_breakup=(PUtils::sampleFlat() < 0.5 ? 0 : 1);
	}
//	cout << beam->GetValue(IS_BREAKUP, &is_breakup) << ":" 
//	     << target->GetValue(IS_BREAKUP, &is_breakup) << endl;
    }

    PParticle participant;
    if (spectator->ID() == 14) {
	participant.SetID(13);
	participant_mass = makeStaticData()->GetParticleMass("n");
    }
    else {
	participant.SetID(14);
	participant_mass = makeStaticData()->GetParticleMass("p");
    }

    Int_t makeloop=kTRUE;
    num_of_realevents++;

    while (makeloop) {
	num_of_sampledevents++;
	exp_w_mean = (Double_t)num_of_realevents/(Double_t)num_of_sampledevents;
	//cout << exp_w_mean << endl;
	//cout << num_of_realevents << ":" << num_of_sampledevents << endl;
	makeloop=kFALSE;
	if (!composite_model || ((composite_model->GetVersionFlag()) & VERSION_WEIGHTING)) {
	    //no consecutive decay model or 
	    //Taken into account in next step
	    Double_t massS, eS, eP, ptot, px, py, pz, t=-1., mdeut;
	    while (t<0.) {
		ptot = SampleFermi(px,py,pz);                  // Fermi momentum
		massS = spectator->M();                        // mass of spectator nucleon
		eS = sqrt(ptot*ptot + massS*massS);            // spectator total energy in deuteron c.m.
		mdeut = makeStaticData()->GetParticleMass("d");
		t = pow(mdeut-massS,2) - 2.*mdeut*(eS-massS);  // off-shell mass**2 of participant
	    }
	    
	    eP = sqrt(ptot*ptot + t);         // participant total energy
	    
	    participant.SetPxPyPzE(-px,-py,-pz,eP);  // initialize participant nucleon	    
	    spectator->SetPxPyPzE(px,py,pz,eS);      // initialize spectator nucleon
	} else {
	    //here we have to fold in consecutive decay model
	    
	    
	    if (beam_breakup) {
		//beam is (breakup-)deuteron
		my_beam = *target;
		my_beam.Boost(-beam->BoostVector());
	    } else {
		my_beam = *beam;
		my_beam.Boost(-target->BoostVector());
	    }
	    if (!mesh) {
		debug_print=0;
		mesh = new PAdaptiveMeshN(0x3, 2, this, 1.);
		mesh->SetRange(0,0.,0.300); //max fermi momentum
		mesh->SetRange(1,-1.,+1.);  //range for cos theta
		mesh->ReCalcYMax();
		mesh->SetThreshold(1.5,0.01);
		mesh->SetMCPoints(1000);
		mesh->Divide(5,2);
		
	    }
	    debug_print=1;
	    mesh->GetRandom();
	    
//	double p = mesh->GetArrayValue(0);
// 	double theta = acos(mesh->GetArrayValue(1));
	    
// 	cout << mesh->GetArrayValue(0) << ":" << mesh->GetArrayValue(1) << endl;
	    
// 	double phi = 2.*TMath::Pi()*PUtils::sampleFlat();
// 	double sth=sin(theta);
// 	px = p*cos(phi)*sth;
// 	py = p*sin(phi)*sth;
// 	pz = p*cos(theta);
//	cout << px << ":" << py << ":" <<pz << endl;
	    
	    Double_t massS, eS, eP, ptot, t=-1., mdeut;
	    
	    ptot = sqrt(px*px+py*py+pz*pz);
	    
	    massS = spectator->M();                        // mass of spectator nucleon
	    
//	cout << ptot << ":" << massS << endl;
	    eS = sqrt(ptot*ptot + massS*massS);            // spectator total energy in deuteron c.m.
	    mdeut = makeStaticData()->GetParticleMass("d");
	    t = pow(mdeut-massS,2) - 2.*mdeut*(eS-massS);  // off-shell mass**2 of participant
	    eP = sqrt(ptot*ptot + t);         // participant total energy	
	    participant.SetPxPyPzE(px,py,pz,eP);  // initialize participant nucleon
	    spectator->SetPxPyPzE(-px,-py,-pz,eS);      // initialize spectator nucleon
//	participant.Print();
//	exit(1);
	}
	
	
	//Up to now we are in the (breakup-)deuteron frame.
	//Let us go into the lab frame first
	if (beam_breakup) {
	    //beam is (breakup-)deuteron
	    participant.Boost(beam->BoostVector()); //go to Lab frame
	    spectator->Boost(beam->BoostVector());  //go to Lab frame
	    if (target->ID() == p1->ID()) { //Identify the scattered nucleon
		*p1=*target;
		*p2=participant;
	    }
	    else  { //Identify the scattered nucleon
		*p2=*target;
		*p1=participant;
	    }
	    
	} else{
	    //target is (breakup-)deuteron
	    //for consistency (target might move) we boost also here
	    participant.Boost(target->BoostVector()); //go to Lab frame
	    spectator->Boost(target->BoostVector());  //go to Lab frame
	    if (beam->ID() == p1->ID()) { //Identify the scattered nucleon
		*p1=*beam;
		*p2=participant;
	    }
	    else  { //Identify the scattered nucleon
		*p2=*beam;
		*p1=participant;
	    }
	} 
	//go into parent frame
	p1->Boost(-parent->BoostVector());
	p2->Boost(-parent->BoostVector());
	spectator->Boost(-parent->BoostVector());
	composite->Reconstruct(); //reset mass after p1 and p2 have been setted

	//now check if we fall into the kinematic limits
	if (composite_model && ((composite_model->GetVersionFlag()) & VERSION_WEIGHTING)) {

	    if (composite_model->GetWeight(composite->M()) <= 0.) {		
		//cout << composite->M() << endl;
		makeloop=kTRUE;
	    }
	}

    }

    //boost scatter back to lab
    p1->Boost(parent->BoostVector());
    p2->Boost(parent->BoostVector());

    spectator->SetW(parent->W());                  // copy parent weight to spectator

    return kTRUE;
}


Bool_t PFermiMomentum:: SampleMomentum(void) {
    return kTRUE;
}


void PFermiMomentum::Print(const Option_t* delme) const {

    PDistribution::Print();
}


