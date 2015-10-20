/////////////////////////////////////////////////////////////////////
//
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "PHadronDecayM1.h"

PHadronDecayM1::PHadronDecayM1()  {
} ;

PHadronDecayM1::PHadronDecayM1(const Char_t *id, const Char_t *de, Int_t key) :
    PHadronDecay(id, de,key) {
    if (is_channel<0)
	Warning("PHadronDecayM1","The model (%s) should be bound to CHANNELS only",de);
  
    //Get particles
    Int_t tid[11];
    tid[0]=10; 
    makeStaticData()->GetDecayModeByKey(key,tid); // retrieve current mode info
    parent_id=makeStaticData()->GetDecayParentByKey(key);

    if (tid[0]!=2) 
	Warning("PHadronDecayM1","(%s):  Only 2 body decay",de);

    //First determine unstable position
    //The hadron with larger width seems to be unstable

    if (makeStaticData()->GetParticleTotalWidth(tid[2]) >
	makeStaticData()->GetParticleTotalWidth(tid[1])) {
	unstable_position=2;
	stable_position=1;
    } else{
	unstable_position=1;
	stable_position=2;
    }
	

    unstable_mass=makeStaticData()->GetParticleMass(tid[unstable_position]); 
    // unstable particle mass pole
    unstable_ml=PData::LMass(tid[unstable_position]);       
    // unstable particle mass threshold

    stable_mass=makeStaticData()->GetParticleMass(tid[stable_position]); 
    unstable_id=tid[unstable_position];
    stable_id  =tid[stable_position];
    scale=1.1;  // good starting value, by trial & error

    //same order then in data base
    id1=tid[1];
    id2=tid[2];
    parent=NULL;
    daughter1=NULL;
    daughter2=NULL;
    model1=NULL;
    model2=NULL;
    didx1=-1;
    didx2=-1;
    didx_unstable=-1;
} ;

PDistribution* PHadronDecayM1::Clone(const char*delme) const {
    return new PHadronDecayM1((const PHadronDecayM1 &)* this);
};

Bool_t PHadronDecayM1::Init(void) {
    //Init function called once for each PChannel

    //looking for parent. This is mandatory
    parent = GetParticle("parent");
    if (!parent) {
	Warning("Init","Parent not found");
	return kFALSE;
    }

    //getting 2 daughters
    //numbering has same order than in data base
    daughter1 = GetParticle(makeStaticData()->GetParticleName(id1));
    daughter2 = GetParticle(makeStaticData()->GetParticleName(id2));

    return kTRUE;

}

Double_t PHadronDecayM1::EvalPar(const Double_t *x, const Double_t *params) {
    if (params) {
	draw_option=(int)params[0];
	didx_option=(int)params[1];
    }
    return Eval(x[0]);
}
 
Double_t PHadronDecayM1::Eval(Double_t x, Double_t y , Double_t z , Double_t t ) const
{
    Double_t res;
    Double_t mass[3];
    Int_t didx[3];
    didx[0]=didx[1]=didx[2]=-1;

    //first determine if parent mass is fixed or has distribution
    PChannelModel *pmodel = makeDynamicData()->GetParticleModel(parent_id);

    mass[unstable_position]=x;
    mass[stable_position]=stable_mass;
    didx[unstable_position]=didx_option;
    
    if (draw_option==0) {

	if (!pmodel)
	    return ((PChannelModel*)this)->GetWeight(mass,didx);
	Double_t w=0;
	Double_t dm=PData::UMass(parent_id)-PData::LMass(parent_id);
	for (mass[0]=PData::LMass(parent_id);mass[0]<(PData::LMass(parent_id)+dm);
	     mass[0]+=dm/100.) { //Convolution over parents distribution

	    w+=((PChannelModel*)this)->GetWeight(mass,didx)*pmodel->GetWeight(mass,(Int_t*)&is_channel);
	}
	return w;
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



int PHadronDecayM1::GetDepth(int i) {
    //Initialize daughter decay modes

    Int_t a1,a2;
    model1 = makeDynamicData()->GetParticleModel(id1);
    model2 = makeDynamicData()->GetParticleModel(id2);

    if (model1) a1=model1->GetDepth(i);
    if (model2) a2=model2->GetDepth(i);

    //after the GetDepth, the threshs are initialized
    makeStaticData()->SetDecayEmin(is_channel, 
				   TMath::Min(
				       makeStaticData()->GetParticleEmin(id1),
				       makeStaticData()->GetParticleEmin(id2)));
    
    return TMath::Max(a1+1,a2+1); 
}


Bool_t PHadronDecayM1::Prepare(void) {
    //Things which might change during the eventloop
    
    didx1=daughter1->GetDecayModeIndex(1);
    didx2=daughter2->GetDecayModeIndex(1);

    //set also unstable id (see sampleM1 why!)
    if (unstable_position==1)
	didx_unstable=didx1;
    else
	didx_unstable=didx2;

    return kTRUE;

}

void PHadronDecayM1::SubPrint(Int_t opt) const {
    //Print sub-models    
    if (model1) {cout << " ";cout << model1->GetDescription();}
    if (model2) {cout << " ";cout << model2->GetDescription();}
}

Double_t PHadronDecayM1::GetWeight(void) {
    Double_t mass[3];
    mass[0]=parent->M();
    mass[1]=daughter1->M();
    mass[2]=daughter2->M();
    return GetWeight(mass); //didx already set in Prepare()

}

Double_t PHadronDecayM1::GetWeight(Double_t *mass, Int_t *didx) {
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

    if (unstable_position==1)
	return BWWeight(unstable_id, mass[0],  mass[1] , mass[2], didx_local1, 0 , didx_local2); 
    return BWWeight(unstable_id, mass[0],  mass[2] , mass[1], didx_local2, 0 , didx_local1); 
}


Bool_t PHadronDecayM1::SampleMass(void) {
    //Mass-sampling wrapper
    Double_t mass[3];
    mass[0]=parent->M();

    
    if (!SampleMass(mass)) return kFALSE;

    daughter1->SetM(mass[1]);
    daughter2->SetM(mass[2]);

    return kTRUE;
};

Bool_t PHadronDecayM1::SampleMass(Double_t *mass, Int_t *didx) {
    //samples the mass of the unstable decay product
    //Input: mass[0] (parent)
    //Output: mass[1] and mass[2] 
    //Order of particles in the array must be the same as defined in data base

    
    if (didx) didx_unstable=didx[unstable_position]; //overwrite what has been set in Prepare()

    Double_t ecm=mass[0];    
    if (sampleM1(ecm)) {
	mass[unstable_position]=unstable_dynamic_mass;
	mass[stable_position]=stable_mass;
	
    }
    else {
	Warning("SampleMass","failed, ecm=%f",ecm);
	return kFALSE;
    }
    if (ecm < (stable_mass+unstable_dynamic_mass)) {
	cout << "PHadronDecayM1::sampleMass: violation of energy" << endl;
	return kFALSE;
    }

    return kTRUE;
};

Bool_t PHadronDecayM1::GetWidth(Double_t mass, Double_t *width, Int_t didx) {

    if (makeStaticData()->GetPWidx(is_channel)==-1) return 0.; // Disabled --> BUGBUG why not static?

    if (!makeStaticData()->GetPWidx(is_channel)) { // Enter below only on the first call

	Info("GetWidth","Called for %s", makeStaticData()-> GetDecayName(is_channel));

	makeDynamicData()->GetParticleDepth(parent_id); // if 1st call will initialize flags

	mmin=makeStaticData()->GetDecayEmin(is_channel);  // mass threshold for the decay
	Double_t w0=makeStaticData()->GetDecayBR(is_channel);      // starting weight
	
	mmax=PData::UMass(parent_id);                // mass ceiling
	Double_t dm=(mmax-mmin)/(maxmesh-3.);          // mass increment
	double mass_threshold, mass_ceiling;
	mass_threshold=PData::LMass(unstable_id);
	mass_ceiling  =mmax-PData::LMass(stable_id);
    
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
	    for (int mc=0;mc<mc_max;mc++) {
		Double_t running_unstable_mass=mass_threshold+((Double_t(mc)/Double_t(mc_max))
							       *(mass_ceiling-mass_threshold));
		Double_t w = makeDynamicData()->GetParticleTotalWeight(running_unstable_mass ,unstable_id);
		if (w>bw_max) bw_max=w;
	    }


	    //Continue with the integration, when a pole mass has been found
	    if (bw_max>0) 
		for (int mc=0;mc<mc_max;mc++) {
		    Double_t temp1 = 0.;

		    //We integrate always over the same binning structure,
		    //this avoids artefacts (IF, 6.6.2007)

		    Double_t running_unstable_mass=mass_threshold+((Double_t(mc)/Double_t(mc_max))
								   *(mass_ceiling-mass_threshold));
		    if (mm>(running_unstable_mass+stable_mass)) {
			//Fold with pcms (phase space factor for BW) 
			temp1=PKinematics::pcms(mm,running_unstable_mass, stable_mass);
			//Fold with mass shape of unstable particle
			Double_t w = makeDynamicData()->GetParticleTotalWeight(running_unstable_mass ,unstable_id);
			temp1*=w;
			if (w>(0.01*bw_max)) {  
			    //Cut-off condition with avoids arteficialy destructed behaviour
			    temp0_norm+=temp1;
			    //Get the Gamma_m_m1_m2			    
			    temp1*=HadronWidth(mm,running_unstable_mass, stable_mass);
			    temp0+=temp1;
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



bool PHadronDecayM1:: sampleM1(const double & ecm) {

    // Mass sampling algorithm in case of hadronic decay to two products
    // one of which is unstable and the other stable
    // Arguments: parent cm energy, 
    // 
    // Data members usage:
    // unstable-product id, unstable-product
    // mass "unstable_mass" (to be returned), stable-product mass

    double unstable_mh = ecm-stable_mass, // unstable particle maximum available mass
	ma, peak, a1, a2, area, b=-1, max1;
    unstable_dynamic_mass=0.;       // reset mass

    if (unstable_mh<=unstable_ml) {
	cout << "PHadronDecayM1:: sampleM1: not enough energy" << endl;
	return kFALSE;
    }

 repeat:                            // re-enter if parameter change is forced
    if (unstable_mh>unstable_mass) {  // mass pole included in the mass range

	peak = BWWeight(unstable_id,scale*ecm,unstable_mass,stable_mass, didx_unstable);

	// probability distribution maximum
	
	ma = 0.5*(unstable_mass+unstable_mh);                // test mass to the right of the pole
	a2 = 0.5*peak*(unstable_mh-ma);           // area in Region II, ma < mass <= mh
	b = peak/(unstable_mh - ma);              // parameters for Region II
    } else {                           // the mass pole is above available energy


	peak = scale*maxBWWeight(unstable_id,ecm,unstable_ml,unstable_mh,max1,stable_mass, didx_unstable);
	// Get max value of mass curve
	ma = unstable_mh;                         // left tail of the distribution function
	a2 = 0.;                         // only Region I, no Region II	
    }
    if (peak==0.) { 

	//if peak not found -> do not sample mass
	unstable_dynamic_mass = 0.;
	return kFALSE;
    }
    a1 = peak*(ma-unstable_ml);                 // area in Region I: unstable_ml < mass <= ma
    area = a1 + a2;                    // test-function total integral
  
    //__________________________________________________________________
    // The rejection method is used for mass sampling:
    double a, f, y, da;                // f(m) is the test function
    // The test function must be integrable, invertible, and bounding the actual
    // distribution function from above over the range of the variable (mass);
    // see e.g. Ref 6. This is by far the most efficient method in this case.
    // Using TF2 and GetRandom of ROOT is several orders of magnitude slower!
  
    do {                               // enter the rejection-method loop
	a = area * PUtils::sampleFlat(); // sample a test area 0 < a(m) < total area
    
	if (a<=a1) {    // Region I:  
	    f = peak;   // test function f(m) is constant over unstable_ml < m <= ma
	    unstable_dynamic_mass = a/f + unstable_ml;  
	    // invert Integral{ f(m) }_unstable_ml^m = a to solve for m
	} else {        // Region II:
	    da = 2.*(area-a);
	    f = sqrt(da*b);                
	    // f(m) is the line from f(ma)=peak to f(unstable_mh)=0 over
	    unstable_dynamic_mass = unstable_mh - da/f;     
	    // ma < m <= unstable_mh; solve Integral{ f(m) }_ma^m = a-a1 for m
	}
    
	y = BWWeight(unstable_id,ecm,unstable_dynamic_mass,stable_mass, didx_unstable);         

    } while (f*PUtils::sampleFlat()>y);// compare with f(m) and accept or reject
  
    if (f<y) {         // Error: the test function fails to bound the distribution fn;
	scale = scale*y/f + 0.1;         // "scale" is too low; adapt (a few events
	// can be done wrong until scale is ok)

	goto repeat;                     // try again
    }
    return kTRUE;
}


double PHadronDecayM1:: BWWeight(const int & i1, const double & m, 
				 const double & m1, const double & m2, int didx_local1, 
				 int i2, int didx_local2) {
    // If i2=0 (default):
    // Breit-Wigner distribution convoluted with a phase-space
    // factor for sampling the mass of an unstable hadron p1
    // in a decay "hadron -> p1 + p2", where p2 is a stable hadron.
    //
    // If i2!=0:
    // "double" Breit-Wigner distribution convoluted with a
    // phase-space factor for the simultaneous sampling of
    // masses of two unstable hadrons p1, p2 in a decay of type
    // "hadron -> p1 + p2"
    // If the correponding didx flags are set (>=0), the BW is taken with 
    // the target decay index
  
    double pCms=PKinematics::pcms(m,m1,m2);

    PChannelModel *i1_model = makeDynamicData()->GetParticleModel(i1);
    if (!i1_model) {
	Warning("BWWeight","No target model i1");
	return 0;
    }
 
    double bw=pCms*i1_model->GetWeight((Double_t*)&m1,&didx_local1);
    
    if (i2) { 

	PChannelModel *i2_model = makeDynamicData()->GetParticleModel(i2);

	if (!i2_model) {
	    Warning("BWWeight","No target model i2");
	    return bw;
	}
	bw*=i2_model->GetWeight((Double_t*)&m2,&didx_local2);
    }
    return bw;
}


double PHadronDecayM1:: maxBWWeight(const int &i1, const double &m, const double &lower,
			 const double &upper, double &max1, const double &m2, int didx_local1, 
				    int i2, int didx_local2) {
    //
    // Find maximum of Weight() using golden rule from Ref.[6] (2nd ed., p. 401)
    //
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

    const double R=0.61803399;
    const double C=1.-R;
    const double tol=1.e-3;
    double ax, bx, cx;
    double f1, f2 ,x0, x1, x2, x3;


    ax = lower;
    bx = C*lower+R*upper;
    cx = upper;
    x0 = ax;
    x3 = cx;
    if ( fabs(cx-bx)>fabs(bx-ax) ) {
	x1 = bx;
	x2 = bx+C*(cx-bx);
    } else {
	x2 = bx;
	x1 = bx-C*(bx-ax);
    }
    f1 = BWWeight(i1, m, x1, m2, didx_local1, i2);
    f2 = BWWeight(i1, m, x2, m2, didx_local1, i2);
    while (fabs(x3-x0)>tol*(fabs(x1)+fabs(x2))) {

	//BUGBUG: No break condition here!
	if (f2>f1) {
	    SHFT3(x0,x1,x2,R*x1+C*x3);
	    SHFT2(f1,f2,BWWeight(i1, m, x2, m2, didx_local1, i2, didx_local2));
	} else {
	    SHFT3(x3,x2,x1,R*x2+C*x0);
	    SHFT2(f2,f1,BWWeight(i1, m, x1, m2, didx_local1, i2, didx_local2));
	}
    }


    if (f1>f2) {
	max1 = x1;
	return f1;
    } else {
	max1 = x2;
	return f2;
    }
}




ClassImp(PHadronDecayM1)
