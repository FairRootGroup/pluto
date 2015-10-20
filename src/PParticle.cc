////////////////////////////////////////////////////////
//  Particle Class implementation file
//
//  A PParticle is a TLorentzVector with id and weight.
//  Static particle properties from the permanent data
//  base in PStaticData are used. Additional external particles
//  and decay modes may be loaded. See also PData.
//
//                    Author:  M.A. Kagarlis
//                    Written: 15.11.98
//                    Revised: 15.12.2000  R. Holzmann
//                    Revised: 22.03.2005  R. Holzmann
//                    Revised: 23.07.2007  IF (new framework)
//
////////////////////////////////////////////////////////

#include "TF1.h"
#include "TMath.h"
#include "PUtils.h"
#include "PParticle.h"
#include "PKinematics.h"
#include "PDynamicData.h"


void PParticle::defaults(void) {
  
    SetVertex(*makeStaticData()->GetBatchValue("_event_vertex_x"),
	      *makeStaticData()->GetBatchValue("_event_vertex_y"),
	      *makeStaticData()->GetBatchValue("_event_vertex_z"),0.);
    pParticle = NULL;
    qParticle1 = NULL; qParticle2 = NULL;
    index = -1;
    parentIndex = -2;
    decayModeIndex = -1;
    destroyDecayModeIndex= 0;
    decayTime = 0.;
    SetSibling(NULL);
    debug="";
    values=NULL;
    make_new_qParticle=1;
    SetSourceId(-1);
    SetSiblingIndex(-1);
    SetDaughterIndex(-1);
    SetParentId(-1);
    SetParentIndex(-1);
    mult=1.;
    spectator=0;
    genwt=invgenwt=1.;

    ResetDaughters();
}

PParticle::PParticle(int id, Double_t T, Double_t w):
    TLorentzVector(0,0,sqrt(T*T+2*T*makeStaticData()->GetParticleMass(id)),
		   T+makeStaticData()->GetParticleMass(id)),
    pid( makeStaticData()->IsParticleValid(id) ), wt( (pid)?w:0.), active(kTRUE) {
    // id, lab kinetic energy (GeV), weight

  defaults();

}


PParticle::PParticle(const char * id, Double_t T, Double_t w):
    TLorentzVector(0,0,sqrt(T*T+2*T*makeStaticData()->GetParticleMass(id)),
		   T+makeStaticData()->GetParticleMass(makeStaticData()->IsParticleValid(id))),
    pid( makeStaticData()->IsParticleValid(id) ), wt( (pid)?w:0.), active(kTRUE) {
    // name, lab kinetic energy (GeV), weight
 defaults();
}


PParticle::PParticle(int id, Double_t px, Double_t py, Double_t pz, Double_t m, Double_t w):
    TLorentzVector(px,py,pz,sqrt(px*px+py*py+pz*pz+
				 ( (m>0.)?m*m:makeStaticData()->GetParticleMass(id)*
				   makeStaticData()->GetParticleMass(id) ) )),
    pid( makeStaticData()->IsParticleValid(id) ), wt( (pid)?w:0.), active(kTRUE) {
    // id, Px, Py, Pz (GeV/c), mass (GeV/c**2) overrides default, weight
 defaults();
}


PParticle::PParticle(const char * id, Double_t px, Double_t py, Double_t pz, Double_t m, Double_t w):
    TLorentzVector(px,py,pz,sqrt(px*px+py*py+pz*pz+
				 ( (m>0.)?m*m:makeStaticData()->GetParticleMass(id)*
				   makeStaticData()->GetParticleMass(id) ) )),
    pid( makeStaticData()->IsParticleValid(id) ), wt( (pid)?w:0.), active(kTRUE) {
    // name, Px, Py, Pz (GeV/c), mass (GeV/c**2) overrides default, weight
 defaults();
}


PParticle::PParticle(int id, const TVector3 & p, Double_t m, Double_t w):
    TLorentzVector( p, sqrt(p.Mag2()+( (m>0.)?m*m:makeStaticData()->GetParticleMass(id)*
				       makeStaticData()->GetParticleMass(id) ))),
    pid(makeStaticData()->IsParticleValid(id) ), wt( (pid)?w:0.), active(kTRUE) {
    // id, 3-momentum vector (GeV/c), mass (GeV/c**2) overrides default, weight
 defaults();
}


PParticle::PParticle(int id, Double_t * pt, Double_t w):
    TLorentzVector( pt ), pid( id ), wt( w ), active(kTRUE) {
    // id, pointer to 4-dim array (Px, Py, Pz, E) (GeV/c, GeV), weight
 defaults();
}

  
PParticle::PParticle(int id, float * pt, Double_t w):
    TLorentzVector( pt ), pid( id ), wt( w ), active(kTRUE) {
    // id, pointer to 4-dim array (Px, Py, Pz, E) (GeV/c, GeV), weight
 defaults();
}


PParticle::PParticle(const PParticle & p):
    TLorentzVector( p.Vect4() ), pid(p.ID() ),
    sourceId( p.GetSourceId() ),
    parentId( p.GetParentId() ),
    parentIndex( p.GetParentIndex() ),
    daughterIndex( p.GetDaughterIndex() ),
    siblingIndex( p.GetSiblingIndex() ),
    decayTime( p.GetProperTime() ), 
    wt( p.W() ), active( p.IsActive() ) {
    // copy constructor
    SetVertex(p.X(),p.Y(),p.Z(),p.T());
    pParticle = NULL;
    qParticle1= NULL;
    qParticle2= NULL;

    make_new_qParticle=p.make_new_qParticle;
    if (p.make_new_qParticle) {
	qParticle1 = (p.qParticle1 ? 
		      new PParticle(p.qParticle1) : NULL); //Copy, because it will destroy in dtor
    } 

    if (p.make_new_qParticle) {
	qParticle2 = (p.qParticle2 ?
		      new PParticle(p.qParticle2) : NULL);
    }

    sParticle  = p.sParticle;
    //      if (p.debug) 
// 	  debug=p.debug;
//       else
    debug="";
    values=NULL;
    destroyDecayModeIndex= p.destroyDecayModeIndex;
    decayModeIndex= p.decayModeIndex;

    if (p.values) values=new PValues(*(p.values));
    mult=p.mult;
    spectator=p.spectator;
    genwt=p.genwt;
    invgenwt=p.invgenwt;
}


PParticle::PParticle(const PParticle * p):
    TLorentzVector( p->Vect4() ), pid( p->ID() ), 
    sourceId( p->GetSourceId() ),
    parentId( p->GetParentId() ),
    parentIndex( p->GetParentIndex() ),
    daughterIndex( p->GetDaughterIndex() ),
    siblingIndex( p->GetSiblingIndex() ),
    decayTime( p->GetProperTime() ), 
    wt( p->W() ), active( p->IsActive() ) {
    // copy constructor
    SetVertex(p->X(),p->Y(),p->Z(),p->T());
    pParticle = NULL;
    qParticle1= NULL;
    qParticle2= NULL;

    make_new_qParticle=p->make_new_qParticle;
    if (p->make_new_qParticle) {
	qParticle1 = (p->qParticle1 ?
		      new PParticle(p->qParticle1): NULL); //Copy, because it will destroy in dtor
    } 

    if (p->make_new_qParticle) {
	qParticle2 = (p->qParticle2 ?
		      new PParticle(p->qParticle2): NULL);
    }

    sParticle  = p->sParticle;
 
    //      if (p->debug) debug=p->debug;
//       else 
    debug="";
    values=NULL;
    destroyDecayModeIndex= p->destroyDecayModeIndex;
    decayModeIndex= p->decayModeIndex;

    if (p->values) values=new PValues(*(p->values));

    mult=p->mult;
    spectator=p->spectator;
    genwt=p->genwt;
    invgenwt=p->invgenwt;
}


void PParticle:: SetKE(Double_t T) {
    // reset by kinetic energy
    Double_t m=M(), th=Theta(), ph=Phi(), 
	ptot=sqrt(T*T+2*T*m), e=T+m, pxy=ptot*sin(th),
	px=pxy*cos(ph), py=pxy*sin(ph), pz=ptot*cos(th);
    SetPxPyPzE(px,py,pz,e);
}

void PParticle:: SetM(Double_t m) {
    // reset mass
    Double_t m1, px=Px(), py=Py(), pz=Pz();
    if (m>0.) m1=m;            // if non-zero reset mass to argument value
    else {                     // else sample from Breit-Wigner if width > 1 MeV
	Double_t g0= makeStaticData()->GetParticleTotalWidth(pid), 
	    m0=makeStaticData()->GetParticleMass(pid);
	if (g0<1.e-3) m1=m0;     // else fix to centroid value
	else do { m1 = PUtils::sampleBW(m0,g0); } while (m1<0.);
    }
    SetXYZM(px,py,pz,m1);
}

void PParticle:: SetMom(Double_t mom) {
    // reset momentum
    Double_t scaling_factor = mom/Rho();
    Double_t px=Px(), py=Py(), pz=Pz();

    px *= scaling_factor;
    py *= scaling_factor;
    pz *= scaling_factor;

    SetXYZM(px,py,pz,makeStaticData()->GetParticleMass(ID()));
}


Double_t PParticle:: Life(Double_t m, int idx) {
    // lifetime in lab frame
    Double_t r, tau=Gamma()*makeDynamicData()->GetParticleLife(pid,m,idx);
    do { r = PUtils::sampleFlat(); } while (r==0.);
    return -tau*log(r);
}    

void PParticle:: SetProperTime() {
    // lifetime in particle's frame  (in mm/c)
 
    Double_t r, tau;

    if (makeStaticData()->GetParticleTotalWidth(pid) > 0) {
	tau = 6.582122e-25 / makeStaticData()->GetParticleTotalWidth(pid);  // hbar units (GeV s)
	do { r = PUtils::sampleFlat(); } while (r==0.);
	decayTime = -tau*log(r)*3.0e11;           // go from sec to mm/c
	return;
    }
    else if (pid==50 || makeStaticData()->IsParticle(pid,"dilepton")) { // dilepton decays instantly
	decayTime = 0.0;
	return;
    }
    else if (pid>1000) { //composites as well
	decayTime = 0.0;
	return;
    }
    else decayTime = 1.e16;
}



Double_t PParticle:: InvariantT(Double_t m3, Double_t m4, Double_t cos_th_cm) {
    // The Mandelstam invariant t in 1 + 2 --> 3 + 4 scattering,
    // where cos_theta_cm must be between -1 and 1.
    // Use the same function for the Mandelstam variable u, with the
    // changes:  m3 <--> m4 & cos_th_cm -> -cos_th_cm

    if (!qParticle1 || !qParticle2 ) {
	Warning("InvariantT","No beam & target found");
	return 0;	
    }

    //Boost into target c.m. frame
    PParticle beam=*qParticle1;
    PParticle target=*qParticle2;
    PParticle parent=*this;
    beam.Boost(-target.BoostVector());
    parent.Boost(-target.BoostVector());

    Double_t ecm=parent.M();

    Double_t p1=(beam.Vect()).Mag();    // the beam momentum
    Double_t m1s=beam.M2();             // the beam mass squared
//  if (pdNNN==2){p1 /= 2.; m1s /= 4.;}// beam is deuteron
    Double_t m2=target.M();             // the target mass (target stationary)
//  if (pdNNN==1) m2 /= 2.;            // target is deuteron
// Quasi-free treated explicetly

    Double_t pcm1=p1*m2/ecm;             // cm momentum of 1 & 2
    Double_t pcm3=PKinematics::pcms(ecm,m3,m4);// cm momentum of 3 & 4
    Double_t ecm1=sqrt(m1s+pcm1*pcm1);   // cm energy of 1
    Double_t ecm3=sqrt(m3*m3+pcm3*pcm3); // cm energy of 3
    Double_t t = m1s - 2.*ecm1*ecm3 + 2.*pcm1*pcm3*cos_th_cm + m3*m3;
    return t;
}


PParticle & PParticle::AddTmp( const PParticle & p) {
    // Like "+" but only for tmp reasons in the evtloop
    // No construction of beam+target objects
    SetPxPyPzE(Px()+p[0],Py()+p[1],Pz()+p[2],E()+p[3]);
    wt*=p.W();
    pid=0;
    active=active&p.IsActive();
    return *this;
}

void PParticle::Reconstruct(void) {
    //reconstruct composite from given beam+target

    if ((qParticle1==NULL) || (qParticle2==NULL)) return;

    SetPxPyPzE(qParticle1->Px()+qParticle2->Px(),
	       qParticle1->Py()+qParticle2->Py(),
	       qParticle1->Pz()+qParticle2->Pz(),
	       qParticle1->E()+qParticle2->E());
}

PParticle & PParticle::operator -= ( const PParticle & p) {
    SetPxPyPzE(Px()-p[0],Py()-p[1],Pz()-p[2],E()-p[3]);
    wt*=p.W();
    pid=0;
    active=active&p.IsActive();
    return *this;
}

PParticle & PParticle::operator += ( const PParticle & p) {
    //Make a composite particle
    //This requires some checks, keeping beam+target
    //For loop calculations use the AddTmp()-method

    if (pid/1000>0) {
	Warning("operator+","No more than two constituent particles");
	return *this;
    }
    int id2=1000*p.ID();
    if (makeStaticData()->IsParticleValid(ID())==0 ||
	makeStaticData()->IsParticleValid(p.ID())==0) {
	printf("PParticle:: cannot add invalid particles\n");
	return *this;
    }

    make_new_qParticle=kTRUE;
    qParticle1 = new PParticle(this); 
    qParticle2 = new PParticle(p);
    Info("operator+","(%s) Keeping beam and target particle for further reference", PRINT_AUTO_ALLOC);
    make_new_qParticle=kTRUE;
    if (makeDataBase()->GetEntryInt("pid", id2+ID())<0) {
	Int_t len=strlen(makeStaticData()->GetParticleName(ID()))+
	    strlen(makeStaticData()->GetParticleName(p.ID()))+6;
	char * delme=new char[len];
	
	sprintf(delme,"%s + %s",makeStaticData()->GetParticleName(ID()),
		makeStaticData()->GetParticleName(p.ID()));
	if (!makeStaticData()->AddParticle(id2+ID(), delme,M()+p.M())) {
	    Warning("operator+","Could not add composite");
	} else {
	    char * delme2=new char[len];
	    //Without spaces
	    sprintf(delme2,"%s+%s",makeStaticData()->GetParticleName(ID()),
		    makeStaticData()->GetParticleName(p.ID()));
	    makeStaticData()->AddAlias(delme,delme2);
	    Info("operator+","(%s) The composite %s has been added", 
		 PRINT_AUTO_ALLOC,makeStaticData()->GetParticleName(id2+ID()));
	}
    }

    SetSibling(NULL);
    SetPxPyPzE(Px()+p[0],Py()+p[1],Pz()+p[2],E()+p[3]);

    wt*=p.W();

    pid+=id2;


    active=active&p.IsActive();
    debug="";
    return *this;
}

void PParticle::Scatter(PParticle *p1, PParticle *p2) {
    //Like the "+" but without making new particles

    if (p1->pid/1000>0 || p2->pid/1000>0) {
	Warning("Scatter","No more than two constituent particles");
	return;
    }

    int id1=p1->ID();
    int id2=p2->ID();
    if (makeStaticData()->IsParticleValid(id1)==0 ||
	makeStaticData()->IsParticleValid(id2)==0) {
	printf("PParticle:: cannot add invalid particles\n");
	return;
    }

    make_new_qParticle=kFALSE;
    qParticle1 = p1; 
    qParticle2 = p2;
    pid =id1 + id2*1000;

    if (makeDataBase()->GetEntryInt("pid", pid)<0) {
	Int_t len=strlen(makeStaticData()->GetParticleName(id1))+
	    strlen(makeStaticData()->GetParticleName(id2))+6;
	char * delme=new char[len];
	
	sprintf(delme,"%s + %s",makeStaticData()->GetParticleName(id1),
		makeStaticData()->GetParticleName(id2));
	if (!makeStaticData()->AddParticle(pid, delme,p1->M()+p2->M())) {
	    Warning("Scatter","Could not add composite");
	} else {
	    char * delme2=new char[len];
	    //Without spaces
	    sprintf(delme2,"%s+%s",makeStaticData()->GetParticleName(id1),
		    makeStaticData()->GetParticleName(id2));
	    makeStaticData()->AddAlias(delme,delme2);
	    Info("Scatter","(%s) The composite %s has been added", 
		 PRINT_AUTO_ALLOC,makeStaticData()->GetParticleName(pid));
	}
    }

    SetSibling(NULL);
    SetPxPyPzE(p1->Px()+p2->Px(),p1->Py()+p2->Py(),p1->Pz()+p2->Pz(),p1->E()+p2->E());

    wt = p1->W() * p2->W();

    active = p1->active & p2->IsActive();
    debug="";
}


PParticle & PParticle::operator = ( const PParticle & p ) {
    // assignment
    
    SetVect4(p.Vect4());
    wt = p.W();
    pid = p.ID();
    parentId = p.GetParentId();
    sourceId = p.GetSourceId();
    parentIndex = p.GetParentIndex();
    daughterIndex= p.GetDaughterIndex();
    siblingIndex = p.GetSiblingIndex();
    decayTime = p.GetProperTime();
    active = p.IsActive();
    SetVertex(p.X(),p.Y(),p.Z(),p.T());
//  if (p.debug) debug=p.debug;
  
    // else 
    //debug="";

    pParticle = NULL;
    qParticle1  = NULL;
    qParticle2  = NULL;
    make_new_qParticle = p.make_new_qParticle;

    if (!qParticle1 && p.make_new_qParticle) {
	qParticle1 = (p.qParticle1 ?
		      new PParticle(p.qParticle1): NULL); //Copy, because it will destroy in dtor
    } //do not make new, but rather overwrite the beam+target particle
    else if (p.qParticle1 && p.make_new_qParticle) *qParticle1=*p.qParticle1;

    if (!qParticle2 && p.make_new_qParticle) {
	qParticle2 = (p.qParticle2 ?
		      new PParticle(p.qParticle2): NULL);
    } //do not make new, but rather overwrite the beam+target particle
    else if (p.qParticle2 && p.make_new_qParticle) *qParticle2=*p.qParticle2;

    sParticle  = p.sParticle;
    destroyDecayModeIndex= p.destroyDecayModeIndex;

    debug="";
    values=NULL;
    return *this;
}

PParticle & PParticle::operator *= (const TRotation & m) {
    // multiplication by rotation matrix
    TLorentzVector v=Vect4();
    v*=m;
    SetVect4(v);
    return *this;
}

PParticle & PParticle::operator *= (const TLorentzRotation & m) {
    // multiplication by Lorentz rotation matrix
    TLorentzVector v=Vect4();
    v*=m;
    SetVect4(v);
    return *this;
}

PParticle & PParticle::Transform(const TRotation & m) {
    // Rotation
    TLorentzVector v=Vect4();
    v.Transform(m);
    SetVect4(v);
    return *this;
}

PParticle & PParticle::Transform(const TLorentzRotation & m) {
    // Lorentz-rotation
    TLorentzVector v=Vect4();
    v.Transform(m);
    SetVect4(v);
    return *this;
}

PParticle* PParticle::Clone(const char*delme) const {
    return new PParticle((const PParticle &)* this);
}

void PParticle::Print(const Option_t* delme) const {
    printf( "  %s (%f,%f,%f;%f) wt = %f, m = %f",
	    (pid/1000==0) ? makeStaticData()->GetParticleName(pid) : "quasi-particle",
	    Px(), Py(), Pz(), E(), wt, M() );
    int id1=pid%1000, id2=pid/1000;
    if (id2==0) {
	printf( " pid = %i\n", id1);
    } else {
	printf( " pid1 = %i, pid2 = %i\n", id1, id2);
    }

    if (values) values->Print();

    printf( "  Vertex = %f %f %f\n", fV.X(), fV.Y(), fV.Z());
    if (delme)
	if (strcmp(delme,"scatter") == 0) {
	    cout << "<scatter>" << endl;
	    if (qParticle1) qParticle1->Print();
	    if (qParticle2) qParticle2->Print();
	    cout << "</scatter>:" << endl;
	}
}





ClassImp(PParticle)


