/////////////////////////////////////////////////////////////////////
//
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PHadronDecayM2.h"


PHadronDecayM2::PHadronDecayM2()  {
} ;

PHadronDecayM2::PHadronDecayM2(const Char_t *id,const  Char_t *de, Int_t key) :
    PHadronDecayM1(id, de,key) {
    if (is_channel<0)
	Warning("PHadronDecayM2","The model (%s) should be bound to CHANNELS only",de);

    
} ;

PDistribution* PHadronDecayM2::Clone(const char*delme) const {
    return new PHadronDecayM2((const PHadronDecayM2 &)* this);
};

Bool_t PHadronDecayM2::Init(void) {
    return PHadronDecayM1::Init();
}

int PHadronDecayM2::GetDepth(int i) {    
    //Initialize daughter decay modes
    return PHadronDecayM1::GetDepth(i);
}


Bool_t PHadronDecayM2::Prepare(void) {
    //Things which might change during the eventloop
    return PHadronDecayM1::Prepare();
}

void PHadronDecayM2::SubPrint(Int_t opt) const {
    //Print sub-models    
    if (model1) {cout << " ";cout << model1->GetDescription();}
    if (model2) {cout << " ";cout << model2->GetDescription();}
}

Bool_t PHadronDecayM2::SampleMass(void) {
    //Mass-sampling wrapper
    Double_t mass[3];
    mass[0]=parent->M();
   
    if (!SampleMass(mass)) return kFALSE;

    daughter1->SetM(mass[1]);
    daughter2->SetM(mass[2]);

    return kTRUE;
};

Bool_t PHadronDecayM2::SampleMass(Double_t *mass, Int_t *didx) {
    //samples the mass of the unstable decay products
    //Input: mass[0] (parent)
    //Output: mass[1] and mass[2] 
    //Order of particles in the array must be the same as defined in data base

    //BUGBUG: Why is didx not used here?
 
    Double_t ecm=mass[0];    
    if (sampleM2(ecm)) {
	mass[1]=dynamic_mass1;
	mass[2]=dynamic_mass2;	
    }
    else {
	return kFALSE;
    }
    if (ecm < (dynamic_mass1+dynamic_mass2)) {
	cout << "PHadronDecayM2::sampleMass: violation of energy" << endl;
	return kFALSE;
    }

    return kTRUE;
};

Bool_t PHadronDecayM2::GetWidth(Double_t mass, Double_t *width, Int_t didx) {

    if (makeStaticData()->GetPWidx(is_channel)==-1) return 0.; // Disabled --> BUGBUG why not static?

    if (!makeStaticData()->GetPWidx(is_channel)) { // Enter below only on the first call

	Info("GetWidth","Called for %s", makeStaticData()-> GetDecayName(is_channel));

	makeDynamicData()->GetParticleDepth(parent_id); // if 1st call will initialize flags

	mmin=makeStaticData()->GetDecayEmin(is_channel);  // mass threshold for the decay
	Double_t w0=makeStaticData()->GetDecayBR(is_channel);      // starting weight
	
	mmax=PData::UMass(parent_id);                // mass ceiling
	Double_t dm=(mmax-mmin)/(maxmesh-3.);          // mass increment
	double mass_threshold1, mass_ceiling1;
	double mass_threshold2, mass_ceiling2;
	
	mass_threshold1=PData::LMass(id1);
	mass_ceiling1  =mmax-PData::LMass(id2);
	mass_threshold2=PData::LMass(id2);
	mass_ceiling2  =mmax-PData::LMass(id1);
    
	mesh = new PMesh(maxmesh-2,"mesh");
	for (int i=0;i<maxmesh-2;++i) {
	    Double_t mm=mmin+i*dm;                     // current invariant mass

	    Double_t temp0 = 0.;
	    Double_t temp0_norm = 0.;

	    //First find out the maximum of the unstable weight,
	    //since very small weights folded with very high phase space
	    //correction seem to diverge
	    //Anyhow, these stuff will not contribute so much to the integral
	    Double_t bw_max=0;
	    for (int mc1=0;mc1<mc_max;mc1++) {
		for (int mc2=0;mc2<mc_max;mc2++) {
		    Double_t running_unstable_mass1=mass_threshold1+((Double_t(mc1)/Double_t(mc_max))
								   *(mass_ceiling1-mass_threshold1));
		    Double_t running_unstable_mass2=mass_threshold2+((Double_t(mc2)/Double_t(mc_max))
								     *(mass_ceiling2-mass_threshold2));
		    Double_t w = makeDynamicData()->GetParticleTotalWeight(running_unstable_mass1 ,id1);
		    w *= makeDynamicData()->GetParticleTotalWeight(running_unstable_mass2 ,id2);
		    if (w>bw_max) bw_max=w;
		}
	    }


	    //Continue with the integration, when a pole mass has been found
	    if (bw_max>0) 
		for (int mc1=0;mc1<mc_max;mc1++) {
		    for (int mc2=0;mc2<mc_max;mc2++) {
			Double_t temp1 = 0.;

			//We integrate always over the same binning structure,
			//this avoids artefacts (IF, 6.6.2007)

			Double_t running_unstable_mass1=mass_threshold1+((Double_t(mc1)/Double_t(mc_max))
								   *(mass_ceiling1-mass_threshold1));
			Double_t running_unstable_mass2=mass_threshold2+((Double_t(mc2)/Double_t(mc_max))
									 *(mass_ceiling2-mass_threshold2));
			if (mm>(running_unstable_mass1+running_unstable_mass2)) {
			    //Fold with pcms (phase space factor for BW) 
			    temp1=PKinematics::pcms(mm,running_unstable_mass1, running_unstable_mass2);
			    //Fold with mass shape all unstable particles
			    Double_t w = makeDynamicData()->GetParticleTotalWeight(running_unstable_mass1 ,id1);
			    w *= makeDynamicData()->GetParticleTotalWeight(running_unstable_mass2 ,id2);			    
			    temp1*=w;
			    if (w>(0.01*bw_max)) {  
				//Cut-off condition with avoids arteficialy destructed behaviour
				temp0_norm+=temp1;
				//Get the Gamma_m_m1_m2			    
				temp1*=HadronWidth(mm,running_unstable_mass1, running_unstable_mass2);
				temp0+=temp1;
			    }
			}
		    }
		}
	    if (temp0_norm>0) //Normalization
		temp0 /= temp0_norm;
	    temp0 *= w0;
	    mesh->SetNode(i,temp0); 

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

Double_t PHadronDecayM2::GetWeight(Double_t *mass, Int_t *didx) {
    //Get the weight of the decay mass[0]->mass[1]+mass[2]
    //Input: mass[0] (parent)
    //mass[1] and mass[2] 
    //Order of particles must be the same as defined in data base

    //If set, use parameters
    int didx_local1=didx1,didx_local2=didx2;
    if (didx) {
	didx_local1=didx[1];
	didx_local2=didx[2];
    }

    return BWWeight(id1, mass[0],  mass[1] , mass[2], didx_local1, id2 , didx_local2); 
}


bool PHadronDecayM2:: sampleM2(const double & ecm) {
    
    // Mass sampling algorithm in case of hadronic decay to two unstable hadrons

    dynamic_mass1=0.; dynamic_mass2=0.;                      // reset masses
    double m0[2]={makeStaticData()->GetParticleMass(id1),
		  makeStaticData()->GetParticleMass(id2)};// mass poles
    double ml[2]={PData::LMass(id1),PData::LMass(id2)};// mass thresholds
    double mh[2]={ecm-ml[1],ecm-ml[0]};    // maximum available masses
    double f, y, m1, m2, mx[2]={TMath::Min(m0[0],mh[0]),TMath::Min(m0[1],mh[1])};
    double maxW1, maxW2;
    
    for(int iter=0;iter<5;iter++) {   // find 2-dim maximum of Weight() function 
	maxW1 = maxBWWeight(id1,ecm,ml[0],mh[0],mx[0],mx[1],didx1,id2,didx2);
	maxW2 = maxBWWeight(id2,ecm,ml[1],mh[1],mx[1],mx[0],didx2,id1,didx1);
    }
    int counter = 0;
 repeat:                           // re-enter if parameter change is forced
    f = scale*TMath::Max(maxW1,maxW2);
    if (f==0.) return kFALSE;
    
    // The rejection method is used to sample the masses m1 and m2.
    // (see e.g. Ref 6). This is by far the most efficient method in this case.
    // Using TF2 and GetRandom of ROOT is several orders of magnitude slower!
    
    do {                              // enter the rejection-method loop
	m1 = ml[0] + PUtils::sampleFlat()*(mh[0]-ml[0]);
	m2 = ml[1] + PUtils::sampleFlat()*(mh[1]-ml[1]);
	y = BWWeight(id1,ecm,m1,m2,didx1,id2,didx2); // y(m1,m2) is the distribution function
	counter++;
    } while ( (PUtils::sampleFlat()>y/f) && counter<10000);
    
    if (f<y) {   // Error: the test function fails to bound the distribution fn;
	scale = scale*y/f + 0.1;
//  printf("%s %f %i\n",Message[14],scale,counter);
	if (counter<10000) goto repeat; // ...and try again
    }
    dynamic_mass1=m1;
    dynamic_mass2=m2;
    
    return kTRUE;
}



ClassImp(PHadronDecayM2)
