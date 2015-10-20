/////////////////////////////////////////////////////////////////////
//
// Pluto build-in genbod
// Copied from PData (originally from Ref. 5)
//
// The particles named with "corr1" and "corr2" can be
// folded with the mass distribution named "/correlation"
// 
//
//                                  Author: Kagarlis
//                                  Reimplemented: IF
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PGenBod.h"
#include "PKinematics.h"


ClassImp(PGenBod)

PGenBod::PGenBod() {
    Fatal("PGenBod()","Wrong constructor");
} ;

PGenBod::PGenBod(const Char_t *id, const Char_t *de, Int_t key) : 
    PChannelModel (id, de,key) {

    parent=NULL;
    primary=NULL;
    n_daughters=0;
    correlator=NULL;
    weight_max=1.0000001;

    for (int i=0;i<MAX_GENBOD_NUM;i++) {
	daughter[i]=NULL;
    }
};

PDistribution* PGenBod::Clone(const char*delme) const {
    return new PGenBod((const PGenBod &)* this);
};

Bool_t PGenBod::Init(void) {
    parent = GetParticle("parent");
    if (!parent) {
	Warning("Init","Parent not found");
	return kFALSE;
    }

    daughter[0] = GetParticle("corr1");
    daughter[1] = GetParticle("corr2");
    
    bool myloop=1;
    for (n_daughters=0;n_daughters<MAX_GENBOD_NUM && myloop;n_daughters++) {
	if (!daughter[n_daughters]) daughter[n_daughters]=GetParticle("daughter");
	if (!daughter[n_daughters]) {
	    myloop=0;
	}
    }
    n_daughters--;

    correlator=
	GetSecondaryModel("correlation");
    //if (correlator) SetVersionFlag(VERSION_NO_PHASE_SPACE);

    return kTRUE;    
};

Double_t PGenBod::GetWeight() {

  if (GetVersionFlag(VERSION_WEIGHTING)) return local_weight;

  return 1;

}


Bool_t PGenBod::SampleMomentum(void) {

    double ecm=parent->M();

    const long double twopi=2.*TMath::Pi();
    int i, j, k, i4, nm1=n_daughters-1, nm2=n_daughters-2, nnm4=3*n_daughters-4, ir=n_daughters-3;
    double tm=0., emmin=0., emmax, bang, etm, wtmax=1., c, s, cb, sb, r,
	esys, beta=-1, gamma=-1, aa, w0, em[n_daughters], emm[n_daughters],
	sm[n_daughters], ems[n_daughters], pd[n_daughters], rno[nnm4], pv[4*n_daughters];
    
    
    for (i=0;i<n_daughters;++i) {
	//	daughter[i]->Print();
	em[i]=daughter[i]->M();
    }

//Pluto-Genbod
    for (i=0;i<n_daughters;++i) {
	ems[i]=(em[i]*em[i]);
	tm+= em[i];
	sm[i]=tm;
    }
    emm[0]=em[0];
    etm=ecm-tm;
    emm[nm1]=ecm;
    emmax=etm+em[0];
    //   cout << "*********" << endl;
    for (i=1;i<n_daughters;++i) {
	emmin += em[i-1];
	emmax += em[i];
	//cout << emmax << ":" << emmin  << ":" << em[i] << endl;
	wtmax *= PKinematics::pcms(emmax,emmin,em[i]);
    }
    if (wtmax==0.) {
	Warning("SampleMomentum","wtmax==0");
	exit(1);
	return kFALSE;    
    }

    Int_t w_counter=0;

 repeat:

    if (correlator) {
	Double_t corr;
	correlator->SampleMass(&corr);
	rno[0] = (corr-sm[1])/etm; 
	//cout << corr << ":" << rno[0] << endl;
    }

 repeat2:
    if (correlator) {
	for (i=1;i<nnm4;++i) rno[i]=PUtils::sampleFlat();
    } else {
	for (i=0;i<nnm4;++i) rno[i]=PUtils::sampleFlat();
    }

    if (nm2>0) {
	if (!correlator)
	    PUtils::dsort(rno,nm2);
	else if (nm2>3)
	    PUtils::dsort(rno+2,nm2-2);
	for (i=1;i<nm1;++i) emm[i]=rno[i-1]*etm+sm[i];
    }
    w0 = 1./wtmax;
    
    for (i=0;i<nm1;++i) {
	pd[i]=PKinematics::pcms(emm[i+1],emm[i],em[i+1]);
	
	if (i>0 || !correlator) 
	    w0  *= pd[i];
	else if (pd[i]==0) w0=0;
    }
    if (w0==0.) {
	if (w_counter < 100) {
	    //make another try...
	    w_counter++;
	    goto repeat2;
	}
	if ((w_counter < 110) && correlator) 
	    goto repeat;
	Warning("SampleMomentum","w0==0");
// 	parent->Print();
// 	for (i=0;i<n_daughters;++i) daughter[i]->Print();
	return kFALSE;
    }
    
    local_weight = 1.;

    if (GetVersionFlag(VERSION_SAMPLING) && !GetVersionFlag(VERSION_NO_PHASE_SPACE)) {
	if (w0 > weight_max) {
	    Warning("SampleMomentum","Weight (%lf) > max (%lf)",w0,weight_max);
	    weight_max = w0*1.1;
	    goto repeat2; 
	}
	if (w0/weight_max<PUtils::sampleFlat()) goto repeat2; // take weight into account (R.H.)
    } else if (GetVersionFlag(VERSION_WEIGHTING)) {
	local_weight = w0;
    } 

    pv[0]=0.;
    pv[1]=pd[0];
    pv[2]=0.;
    for (i=1;i<n_daughters;++i) {
	i4=4*i;
	pv[i4]=0.;
	pv[i4+1]=-pd[i-1];
	pv[i4+2]=0.;
	++ir;
	bang=twopi*rno[ir];
	cb=cos(bang);          // cos(phi)
	sb=sin(bang);          // sin(phi)
	++ir;
	r=rno[ir];
	c=2.*r-1.;             // here comes the scattering angle sampling 

	s=sqrt(1.-c*c);        // sin(theta)
	if (i!=nm1) {
	    esys=sqrt(pd[i]*pd[i]+emm[i]*emm[i]);
	    beta=pd[i]/esys;
	    gamma=esys/emm[i];
	}
	for (j=0;j<=i;++j) {
	    //for (j=1;j<=i;++j) {
	    int jj=4*j, jj1=jj+1, jj2=jj1+1, jj3=jj2+1;
	    aa=( pv[jj ] *  pv[jj ] +
		 pv[jj1] *  pv[jj1] +
		 pv[jj2] *  pv[jj2] );
	    pv[jj3]=sqrt(aa + ems[j]);
	    PKinematics::rotes(c,s,cb,sb,pv+jj);
	    if (i!=nm1) pv[jj1] = gamma*( pv[jj1] + beta * pv[jj3] );
	}
  }
    
    for (i=0;i<n_daughters;++i) {
	k=i+1;
	i4=4*i;
	
	// ______________________________________________________________
	// In the original GENBOD the beam was oriented along the y axis.
	// This is irrelevant so long as theta and phi are isotropic.
	// Since we have introduced anisotropic polar angles,
	// p_y must be interchanged with p_z in order that z be
	// the beam axis. This is equivalent to letting the azimuthal
	// angle phi -> 360-phi. For consistency in case phi anisotropy is
	// also introduced in the future, the sign of p_x has been flipped.
	//                                                     MAK 7.10.99
	//	cout << "set:" << -pv[i4]  <<":"<<  pv[i4+2]  <<":"<<  pv[i4+1]   <<":"<<   pv[i4+3] << endl;
	daughter[i]->SetPxPyPzE(-pv[i4],pv[i4+2],pv[i4+1],pv[i4+3]);
	
	// ______________________________________________________________

	
    }

    return kTRUE;

};

void PGenBod::Print(const Option_t* delme) const {
    BasePrint();
}


