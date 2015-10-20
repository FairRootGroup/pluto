/////////////////////////////////////////////////////////////////////
// 
// Dalitz decay following Ref. 12,15
// 
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PDalitzDecay.h"


ClassImp(PDalitzDecay)

PDalitzDecay::PDalitzDecay()  {
} ;

PDalitzDecay::PDalitzDecay(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de,key) {

    if (is_channel<0)
	Warning("PDalitzDecay","This model should be bound to CHANNELS only");

    //set all needed parameters in advance
    alpha=1./137.0359895;
    mass_pi0 =makeStaticData()->GetParticleMass("pi0");
    mass_n   =makeStaticData()->GetParticleMass("n");
    mass_eta =makeStaticData()->GetParticleMass("eta");
    mass_p   =makeStaticData()->GetParticleMass("p");
    p1   =4.*alpha/(3.*TMath::Pi());
    p2   =0.5*p1;
    p3   =0.5*p2*alpha;
    //--> from PData

    Int_t tid[11];
    tid[0]=10; 
    makeStaticData()->GetDecayModeByKey(primary_key,tid); // retrieve current mode info
    parent_id=makeStaticData()->GetDecayParentByKey(primary_key);

    if (tid[0]!=2) 
	Warning("PDalitzDecay","Only 2 body decay");


    //ee or mumu?
    if (makeStaticData()->IsParticle(tid[1],"dilepton")||
	makeStaticData()->IsParticle(tid[2],"dilepton")) sw=0;
    else if (makeStaticData()->IsParticle(tid[1],"dimuon")||
	     makeStaticData()->IsParticle(tid[2],"dimuon")) sw=1;
    else cout << PRINT_WARNING << de <<
	     ": No dilepton/dimuon" << endl;

    if (makeStaticData()->IsParticle(tid[1],"dilepton")||
	makeStaticData()->IsParticle(tid[1],"dimuon")) {
	dilepton_position=1;
	dilepton_pid=tid[1];
	others_position=2;
	other_pid=tid[2];
    }
    else {
	dilepton_position=2;
	dilepton_pid=tid[2];
	others_position=1;
	other_pid=tid[1];
    }

    if (sw==0) { // e+e-
	mass_e=makeStaticData()->GetParticleMass("e-"); 
	mass_ee=2.*mass_e;
    } 
    else { // mu+mu-
	mass_e=makeStaticData()->GetParticleMass("mu-"); 
	mass_ee=2.*mass_e;
    } //--> from PData
    

 
    mass_x=0.;
    if (makeStaticData()->IsParticle(parent_id,"pi0")||
	makeStaticData()->IsParticle(parent_id,"eta")||
	makeStaticData()->IsParticle(parent_id,"eta'")||
	makeStaticData()->IsParticle(parent_id,"J/Psi")) flag=1;
    else if (makeStaticData()->IsParticle(parent_id,"w")) {
	flag=2;
	mass_x=mass_pi0;
    } 
    else if (makeStaticData()->IsParticle(parent_id,"D0") || 
	     makeStaticData()->IsParticle(parent_id,"NS110")) {
	flag=3;
	mass_x=mass_n;
    } 
    else if (makeStaticData()->IsParticle(parent_id,"D+") || 
	     makeStaticData()->IsParticle(parent_id,"NS11+")) {
	flag=3;
	mass_x=mass_p;
    } 
    else if (makeStaticData()->IsParticle(parent_id,"pn")) { 
	flag=3;
	mass_x=mass_p+mass_n;
    } 
    else if (makeStaticData()->IsParticle(parent_id,"phi")) {
	flag=2;
	mass_x=mass_eta; //BUGBUG-> Are we sure about the mass_x? Needs additional checks
    } else cout << PRINT_WARNING << de <<
	     ": Mode not recognized" << endl;

    //-->used for FDalitz
    pi0        = makeStaticData()->GetParticleID("pi0");
    eta        = makeStaticData()->GetParticleID("eta"); 
    eta_prime  = makeStaticData()->GetParticleID("eta'");
    w          = makeStaticData()->GetParticleID("w"); 
    Delta_0    = makeStaticData()->GetParticleID("D0"); 
    Delta_plus = makeStaticData()->GetParticleID("D+"); 
    phi        = makeStaticData()->GetParticleID("phi");
    S11_0      = makeStaticData()->GetParticleID("NS110"); 
    S11_plus   = makeStaticData()->GetParticleID("NS11+");

    useQED=0;
    flatMD=0;

    if (sw==0) ml=makeStaticData()->GetParticleMass("dilepton");    // dilepton mass threshold
    else ml=makeStaticData()->GetParticleMass("dimuon");

    mass_parent=makeStaticData()->GetParticleMass(parent_id); 
    draw_parent_mass=-1.; //only for draw option

    rejection_flag = (makeStaticData()->IsParticle(parent_id,"eta'")) ? 1 :   // eta'
	( (makeStaticData()->IsParticle(parent_id,"w")) ? 2 :     // w with mass > pole
	  0 );                           // all the other cases
    photon_br = PhotonBR(parent_id);
    integral = NULL;

    formfactor_model=NULL;

} ;

PDistribution* PDalitzDecay::Clone(const char*delme) const {
    return new PDalitzDecay((const PDalitzDecay &)* this);
};

Bool_t PDalitzDecay::Init(void) {
    //Init function called once for each PChannel
    
    //looking for parent. This is mandatory
    parent = GetParticle("parent");
    if (!parent) {
	Warning("Init","Parent not found");
	return kFALSE;
    }

    PParticle * daughter1 = GetParticle("daughter");
    PParticle * daughter2 = GetParticle("daughter");
    if ((daughter1->ID() == makeStaticData()->GetParticleID("dilepton")) ||
	(daughter1->ID() == makeStaticData()->GetParticleID("dimuon"))) {
	other   =daughter2;
 	dilepton=daughter1;
    }
    else {
	other   =daughter1;
 	dilepton=daughter2;
    }
    
    return kTRUE;
}

Bool_t PDalitzDecay::FreezeOut(void) {
    formfactor_model=
	GetSecondaryModel("formfactor");
    return kTRUE;
}

int PDalitzDecay::GetDepth(int i) {
    if (i) return -1; //no Hdepth
    return 0;
}

Bool_t PDalitzDecay::SampleMass(void) {
    //Mass-sampling wrapper
    Double_t masses[3]; //parent + daughter masses
    //It is important that the orders are the same as in data base
    masses[0]=parent->M();
    masses[others_position]=other->M();
    masses[dilepton_position]=dilepton->M();

    SampleMass(&masses[0]); //Call Coupled-Channel wrapper

    other->SetM(masses[others_position]);
    dilepton->SetM(masses[dilepton_position]);
    return kTRUE;
};

Bool_t PDalitzDecay::SampleMass(Double_t *mass, Int_t *didx) {
    sampleMD(mass[0], parent_id, 
	     mass[dilepton_position], mass[others_position]);
    return kTRUE;
};

Double_t PDalitzDecay::GetWeight(void) {
    //Wrapper
    Double_t masses[3]; //parent + daughter masses
    //It is important that the orders are the same as in data base
    masses[0]=parent->M();
    masses[others_position]=other->M();
    masses[dilepton_position]=dilepton->M();

#if 0
    //the following lines are for the normalization of the weight integral
    //as a function of the parents mass
    //otherwise we get artefacts when using the weight option of Pluto
    if (!integral) {
	Int_t maxmesh=200;
	Double_t mmin=PData::LMass(parent->ID());   // mass threshold
	Double_t mmax=PData::UMass(parent->ID());   // mass ceiling
	
	Double_t dm=(mmax-mmin)/(maxmesh-1.); // mass increment for the mesh
	integral = new PMesh(maxmesh,"mesh2");
	integral->SetMin(mmin);
	integral->SetMax(mmax);
	for (int j=0;j<maxmesh;++j) {  // loop over the mesh points [0-200]
	    Double_t local_int=0;
	    Double_t mm=mmin+j*dm;            // update mass on the mesh
	    for (int i=0;i<1000;i++) {
		sampleMD(mm, parent_id, 
			 masses[dilepton_position], masses[others_position]);
		local_int+=dGdM(parent_id,masses[dilepton_position],mm);

	    }
	    integral->SetNode(j,local_int/1000.);
	    //cout << local_int << endl;
	}
	masses[0]=parent->M();
	masses[others_position]=other->M();
	masses[dilepton_position]=dilepton->M();
    }
#endif
    //cout << "INT:"<< integral->GetLinearIP(masses[0]) << endl;
    //return GetWeight(&masses[0])/integral->GetLinearIP(masses[0]); //Call Coupled-Channel wrapper
    //return GetWeight(&masses[0]);

    //Pluto weight should be in units of probability
     return (GetWeight(&masses[0]) / 
	     makeDynamicData()->GetParticleTotalWidth( parent->M(),parent_id));
//    return GetWeight(&masses[0]) / 
//	makeStaticData()->GetParticleTotalWidth(parent_id);
}

Double_t PDalitzDecay::GetWeight(Double_t *mass, Int_t *didx) {
    Double_t dgdm = dGdM(parent_id,mass[dilepton_position],mass[0]);
    //    Double_t dgdm = FDalitz(parent_id,mass[dilepton_position],mass[0]);

    return dgdm;
}


Bool_t PDalitzDecay::GetWidth(Double_t mass, Double_t *width, Int_t didx) {
    if (makeStaticData()->GetPWidx(is_channel)==-1) return 0.; // Disabled --> BUGBUG why not static?

    if (!makeStaticData()->GetPWidx(is_channel) || width_init == 0) { // Enter below only on the first call
#ifdef INGO_DEBUG
	Info("GetWidth","Creating mesh in %s",makeStaticData()-> GetDecayName(is_channel));
#endif
	width_init++;
	makeDynamicData()->GetParticleDepth(parent_id); // if 1st call will initialize flags

	mmin=makeStaticData()->GetDecayEmin(is_channel);  // mass threshold for the decay
	Double_t w0=1.;      // starting weight
	
	if (PData::IsMDalitz(is_channel)) 
	    w0=photon_br;// if meson Dalitz decay...
	// ...then scale by BR of the corresponding real-photon decay (from PBR[] table)
	mmax=PData::UMass(parent_id);                // mass ceiling
	Double_t dm=(mmax-mmin)/(maxmesh-3.);          // mass increment
	double mass_threshold, mass_ceiling;
	mass_threshold=PData::LMass(dilepton_pid);
	mass_ceiling  =PData::LMass(other_pid);
    
	mesh = new PMesh(maxmesh-2,"mesh");
	for (int i=0;i<maxmesh-2;++i) {
	    Double_t mm=mmin+i*dm;   // current invariant mass of the parent resonance

	    Double_t temp0 = 0.;
	    mc_max = 10000;
	    for (int mc=0;mc<mc_max;mc++) 

		temp0+=dGdM(parent_id, mass_threshold+((Double_t(mc)/Double_t(mc_max))
						       *(mm-mass_threshold-mass_ceiling)), mm);
	    temp0 /= mc_max;

	    temp0 *= w0;
	    temp0 *= (mm-mass_ceiling);
	    mesh->SetNode(i,temp0); 
//	    cout << "inv_mass:" << mm <<":"<< temp0 << endl;
	}

	mesh->SetMin(mmin);                 // store mass threshold
	mesh->SetMax(mmax);                 // store mass ceiling
	makeStaticData()->SetPWidx(is_channel,1);  //Enable channel
    } //END first call

    if (mesh)  {
	*width=mesh->GetLinearIP(mass);
	return kTRUE;
    }
    return kFALSE;
}



Double_t PDalitzDecay::EvalPar(const Double_t *x, const Double_t *params) {
    return Eval(x[0]);
}
 
Double_t PDalitzDecay::Eval(Double_t x, Double_t y , Double_t z , Double_t t ) const
{
    //Draw the Dalitz FF as a function of the dilepton mass
    Double_t res;
    Double_t mass[3];

   
    mass[0]=makeStaticData()->GetParticleMass(parent_id); 
    if (draw_parent_mass>0) mass[0]=draw_parent_mass;
    mass[dilepton_position]=x;
    mass[others_position]=mass_x;

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

double PDalitzDecay:: dGdM(const int& id, const double& m, const double& ecm) {
    // Dalitz dilepton dG/dM, the dilepton mass distribution function (Refs. 12, 8)
    // arguments: id=parent hadron, m=dilepton mass, ecm=invariant mass
    // Note: the returned value in the case of Delta0 or Delta+ is the decay width,
    //       but it is the branching ratio for meson Dalitz decays. In the latter
    //       case the real-photon mode static branching ratio scales dG/dM (Ref. 8)
    //       See Ref. 4 eqs. 3.30, 3.31, 3.36, 3.37 for details.
    //
    // (copied and adapted from PData, IF)


    if (m<mass_ee||m>ecm-mass_x) return 0.;
    
    double f1=mass_e/m, f2=f1*f1, a=1.-4.*f2,
	b=1.+2.*f2, y = sqrt(a)*b/m, c, z=-1;

    switch (flag) {
	case 1:                 // pseudoscalar meson
	    c = m/ecm;
	    z = 1.-c*c;
	    z *= p1*z*z;          // Ref 12 Eq. 3.36
		
	    break;
	case 2:                 // vector meson (w)
	    c = ecm*ecm;
	    c = m/(c-mass_x*mass_x);
	    z = 1.+m*c;
	    z *= z;
	    c *= 2.*ecm;
	    z -= c*c;
	    z *= p2*sqrt(z);      // Ref 12 Eq. 3.37
	    break;
	case 3:                 // Delta resonance, N1535 or pn bremsstrahlung
	    z = p3*PKinematics::pcms(ecm,mass_x,m)/(ecm*ecm);
//	    z *= makeStaticData()->GetParticleTotalWidth(id)/0.12; // fix an asymmetry in the treatment of
	    // leptonic vs. hadronic decays
	    // (compare treatment of width in HWidth())
	    break;
    }

    return y*z*FDalitz(id,m,ecm);
}

double PDalitzDecay:: FDalitz(const int & id, const double & m, double ecm) {
    // Matrix element squared for Dalitz decay (Ref 4).
    // Form factor for pseudoscalars

    double f=0., m2=m*m;

    if (id==Delta_0 || id==S11_0 || id==Delta_plus || id==S11_plus) {
	// For Delta, N* we need the matrix element in addition
	Double_t mn=0.;
	if (id==Delta_0 || id==S11_0) 
	    mn=mass_n;
	else if (id==Delta_plus || id==S11_plus) 
	    mn=mass_p;
	double mm=ecm+mn;
	double ff=2.7;
	if (formfactor_model)
	    ff=formfactor_model->GetWeight((Double_t*)&m);
	f = ff*mm/(2.*mn*(mm*mm-m2));
	double md=ecm, md2=md*md, md3=md2*md, md4=md3*md,
	    mn2=mn*mn, mn4=mn2*mn2, m4=m2*m2,
//	sum = (mn-m-md)*(mn+m-md)*         // Eq. 3.54 of Ref 4
//	    sum = (1./9.)*((md-mn)*(md-mn)-m2)*   //from PRC 58, 447
	    sum = ((md-mn)*(md-mn)-m2)*
	    (7.*md4 + 14.*md2*m2 + 3.*m4 + 8.*md3*mn + 2.*md2*mn2 + 
	     6.*m2*mn2 + 3.*mn4);
	return sum*f*f;
    }

    if (useQED) return 1.;   // QED form factor 
    //Self-defined form factor
    if (formfactor_model) {
	double ff=formfactor_model->GetWeight((Double_t*)&m);
	return ff;
    }

    if (id==pi0) {
	f=(1.+5.5*m2);
	return f*f;
    } else if (id==eta) {
	f=1.-m2/0.5184;
	return 1./(f*f);
    } else if (id==eta_prime) 
	return 0.33939776/( (0.5776-m2)*(0.5776-m2) + 0.005776);
    else if (id==w)
	return 0.17918225/( (0.4225-m2)*(0.4225-m2) + 0.000676);
    else if (id==phi) {
	f=1.-m2/(0.77*0.77);
	return 1./(f*f); 
    } else if (id==makeStaticData()->GetParticleID("J/Psi")) 
	return 1.; 
    
    return 0.;

}


void PDalitzDecay:: sampleMD(const double & ecm, const int & id, 
		      double & m, const double & m1) {
    // Mass sampling algorithm in case of Dalitz decay to a dilepton 
    // and the a second particle whose mass is taken to be fixed.
    // Arguments: parent cm energy, parent id, dilepton mass 
    //            (to be returned), other product mass (fixed)
    //
    // sampling does not work for eta->gmumu  ===> adjust the bounding function
    // or find better method! 


    if (flatMD) {    // for test purposes sample a flat mass distribution
	m = ml + (ecm-ml-m1)*PUtils::sampleFlat();
	return;
    }

    //because there is the eta->gmumu bug and I have no idea how
    //the Kagarlis sampling work, I try the ROOT GetRandom()
    //Ingo, 21.1.2009

    if (id == eta_prime) {
	m = GetRandom();
	//	cout << GetIdentifier() << ":" << m << endl;
	return;
    }

    int stored_rejection_flag= rejection_flag;

    if (rejection_flag == 2) //omega
	if (ecm<=mass_parent)
	    // w with mass < pole: go back to flag 0
	    rejection_flag = 0;


    double mx=ecm-m1;                  // dilepton mass ceiling
    double p1=dGdM(id,ml,ecm);      // absolute peak for m=2me
    double b=ml, dw=ecm-mass_parent, area;      // rejection-method parameters
    double p2=-1, mi, a, y, bw=-1, a1=-1, a2, f;


    if (rejection_flag==0 && sw==1) p1=dGdM(id,ml+0.1,ecm);  // for mumu decay add a bit to get peak

    if (rejection_flag==1) b /= ecm;             // adjust parameter for eta'

    if (rejection_flag==2&&ecm>mass_parent) {  // parameter area for w with ecm>pole mass
	mi = mass_parent - m1;       // cusp positioned here if omega with mass > pole
	p2 = (ecm+ml-mi)*dGdM(id,mi,ecm);// value of the local peak
	a1 = (p1+p2)*log(mi/ml);         // area up to m = m_w0 - m_pi0
	bw = p1/mi + p2/(dw+ml);
	a2 = 0.5*bw*dw;                  // area past mass_parent, containing the cusp
	area = a1 + a2;                  // total area under the test function
    } else {                             // area for all the other cases
	area = log(mx/ml);
	if (rejection_flag==1) area *= 2.;         // adjust if eta'
    }
  
    //__________________________________________________________________
    // The rejection method is used for mass sampling:
    // The test function must be integrable, invertible, and bounding the actual
    // distribution function from above over the range of the variable (mass);
    // see e.g. Ref 6.
  
    do {                               // enter the rejection-method loop

	a = area * PUtils::sampleFlat(); // sample a test area 0 < a(m) < total area

	if (rejection_flag!=2) {
	    m = b * exp(a);
	    if (rejection_flag==1) m*=(ecm+ml)/(1.+m);
	    f=p1/m;
	    if (rejection_flag==1) f+=p1/(ecm+ml-m);
	} else if (a>a1) {
	    m = 2.*mx*(area-a)/bw;
	    m = mx - sqrt(m);
	    f=-bw*(m-mx)/dw;
	} else {
	    m = ml*exp(a/(p1+p2));
	    f=p1/m+p2/(mx+ml-m);
	}

	y = dGdM(id,m,ecm);          // y(m) is the actual distribution function
    
    } while (f*PUtils::sampleFlat()>y);// compare with f(m) and accept or reject
  

    rejection_flag = stored_rejection_flag;
    if (f<y) m=ml; //method fails, return minimum mass
}


double PDalitzDecay:: PhotonBR(const int & id) {
  // Looks up (from PBR[]) the static BR for a meson with a Dalitz decay
  // channel X + dilepton to the corresponding X + (real photon) channel.
  // (Look at the comments in dG/dM to see why this is needed.)
    
    int i, np, i1, i2;
    
    int x = 0;   // this is the 2nd product (non dilepton) id
    if (makeStaticData()->IsParticle(id,"pi0")||
	makeStaticData()->IsParticle(id,"eta")||
	makeStaticData()->IsParticle(id,"eta'")) 
	x = makeStaticData()->GetParticleID("g");
    else if (makeStaticData()->IsParticle(id,"D0") || 
	     makeStaticData()->IsParticle(id,"Lambda") || 
	     makeStaticData()->IsParticle(id,"NS110") || 
	     makeStaticData()->IsParticle(id,"ND130")) 
	x = makeStaticData()->GetParticleID("n");
    else if (makeStaticData()->IsParticle(id,"D+")   ||
	     makeStaticData()->IsParticle(id,"NP11+") ||
	     makeStaticData()->IsParticle(id,"ND13+") ||
	     makeStaticData()->IsParticle(id,"NS11+")) 
	x = makeStaticData()->GetParticleID("p");
    else if (makeStaticData()->IsParticle(id,"w")||
	     makeStaticData()->IsParticle(id,"phi")) 
	x = makeStaticData()->GetParticleID("pi0");
    Int_t tid[11];

    Int_t *didx;
    Int_t listkey=-1;
    int key=makeDataBase()->GetEntryInt("pid", id);
    while (makeDataBase()->MakeListIterator
	   (key, "pnmodes", "link", &listkey)) {	
	makeDataBase()->GetParamInt 
	    (listkey, "didx" , &didx); //small workaround -> should work on keys
	tid[0]=10;
	makeStaticData()->GetDecayModeByKey(listkey,tid); // retrieve current mode info
	int pos=makeStaticData()->GetDecayIdxByKey(listkey);
	
	for (i=0;i<tid[0];++i) {            // loop over the decay modes
	    np=tid[0];                    // number of products in current mode
	    if (np==2) {                   // check 2-product decays only
		i1=tid[1];                  // pids of products
		i2=tid[2];
		if ((makeStaticData()->IsParticle(i1,"g")
		     &&i2==x)||    // good decay mode
		    (makeStaticData()->IsParticle(i2,"g")
		     &&i1==x)) return makeStaticData()->GetDecayBR(pos);
	    };
	}
    }
    
    return 0.;                       // no such mode
}


void PDalitzDecay::Print(const Option_t* delme) const {  
    //Debug info
    BasePrint();
    if (useQED) cout << "Uses QED form factor" << endl;
    else cout << "Uses VMD form factor" << endl;
}
