/////////////////////////////////////////////////////////////////////
//
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PEEDirectDecay.h"


ClassImp(PEEDirectDecay)

PEEDirectDecay::PEEDirectDecay()  {
} ;

PEEDirectDecay::PEEDirectDecay(const Char_t *id, const Char_t *de, Int_t key) :PChannelModel(id, de,key) {
    if (is_channel<0)
	Warning("PEEDirectDecay","The model (%s) should be bound to CHANNELS only",de);

    //set all needed parameters in advance

    Int_t tid[11];
    tid[0]=10; 
    makeStaticData()->GetDecayModeByKey(key,tid); // retrieve current mode info
    parent_id=makeStaticData()->GetDecayParentByKey(key);
    
    if (tid[0]!=2) 
	Warning("PEEDirectDecay","(%s) Only 2 body decay",de);

    if (makeStaticData()->IsParticle(parent_id,"rho0"))       
	cv=3.079e-6;   
    else if (makeStaticData()->IsParticle(parent_id,"w"))     
	cv=0.287e-6;
    else if (makeStaticData()->IsParticle(parent_id,"phi"))   
	cv=1.450e-6;
    else if (makeStaticData()->IsParticle(parent_id,"J/Psi")) 
	cv=1.560e-4;
    // eta-> ee: no vector meson 
    else {
	Warning("PEEDirectDecay","(%s) Undetermined direct decay",de);
    }

    //ee or mumu?
    if (makeStaticData()->IsParticle(tid[1],"e+") ||
	makeStaticData()->IsParticle(tid[2],"e+")) mlep=makeStaticData()->GetParticleMass("e-");
    else if (makeStaticData()->IsParticle(tid[1],"mu+")||
	     makeStaticData()->IsParticle(tid[2],"mu+")) mlep=makeStaticData()->GetParticleMass("mu-");

    else Warning("PEEDirectDecay","(%s) No dilepton/dimuon",de);

    use_pi_cutoff=0;
    use_hadronic_ps=0;
};

PDistribution* PEEDirectDecay::Clone(const char*delme) const {
    return new PEEDirectDecay((const PEEDirectDecay &)* this);
};

Bool_t PEEDirectDecay::Init(void) {
    //Init function called once for each PChannel

    //looking for parent. This is mandatory
    parent = GetParticle("parent");
    if (!parent) {
	cout << "PEEDirectDecay::init: parent not found" << endl;
	return kFALSE;
    }

    //getting leptons
    e1 = GetParticle("daughter");
    e2 = GetParticle("daughter");
    return kTRUE;
}

int PEEDirectDecay::GetDepth(int i) {
    makeStaticData()->SetDecayEmin(is_channel, 2*mlep);
    if (i) return -1; //no Hdepth
    return 0;
}

Bool_t PEEDirectDecay::SampleMass(void) {
    //Mass-sampling wrapper
    e1->SetM(mlep);
    e2->SetM(mlep);
    return kTRUE;
};

Bool_t PEEDirectDecay::SampleMass(Double_t *mass, Int_t *didx) {
    //Set to fixed decay products here
    mass[1]=mlep;
    mass[2]=mlep;
    return kTRUE;
};

#if 1
Bool_t PEEDirectDecay::GetWidth(Double_t mass, Double_t *width, Int_t didx) {

    if (makeStaticData()->GetPWidx(is_channel)==-1) {
	*width = makeStaticData()->GetDecayPartialWidth(is_channel);
	return kFALSE;
    }

    if (!makeStaticData()->GetPWidx(is_channel) || width_init == 0) { // Enter below only on the first call
#ifdef INGO_DEBUG
	Info("GetWidth","Creating mesh in %s",makeStaticData()-> GetDecayName(is_channel));
#endif
	width_init++;
	makeDynamicData()->GetParticleDepth(parent_id); // if 1st call will initialize flags

	mmin=makeStaticData()->GetDecayEmin(is_channel);  // mass threshold for the decay
	
	mmax=PData::UMass(parent_id);                // mass ceiling
	Double_t dm=(mmax-mmin)/(maxmesh-3.);          // mass increment
    
	mesh = new PMesh(maxmesh-2,"mesh");
	for (int i=0;i<maxmesh-2;++i) {
	    Double_t mm=mmin+i*dm;   // current invariant mass of the parent resonance

	    Double_t temp0 = 0.;

	    //--------------------
	    double me=2.*mlep*mlep;

	    double m2=mm*mm, m3=m2*mm, mem=me/m2;
	    double wid = 0;
	    if ((1.-2.*mem)>0)
		wid = cv*sqrt(1.-2.*mem)*(1.+mem)/m3;
	    
	    //pipi cutoff
	    double pi2phase = PKinematics::pcms(mm,makeStaticData()->GetParticleMass("pi0"),
						makeStaticData()->GetParticleMass("pi0"));
	    
	    if ((wid>0) && (use_pi_cutoff)) {
		temp0= wid*pow(pi2phase,3./2.) ;
	    } else if ((wid>0) && (use_hadronic_ps)) {
		Double_t pole_width = 
		    makeDynamicData()
		    ->GetParticleTotalWidthSum(makeStaticData()->GetParticleMass(parent_id),parent_id,1);
		Double_t mass_width = 
		    makeDynamicData()->GetParticleTotalWidthSum(mm,parent_id,1);
		
		temp0= wid*(mass_width/pole_width);
	    } else if (wid>0) {
		temp0= wid;
	    } else temp0=0;
	    
	    //--------------------

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
#endif


#if 0
Bool_t PEEDirectDecay::GetWidth(Double_t mass, Double_t *width, Int_t didx) {
    
    if (makeStaticData()->GetPWidx(is_channel)==-1) {
	*width = makeStaticData()->GetDecayPartialWidth(is_channel);
	return kFALSE;
    }
    
    double me=2.*mlep*mlep;

    double m2=mass*mass, m3=m2*mass, mem=me/m2;
    double wid = 0;
    if ((1.-2.*mem)>0)
	wid = cv*sqrt(1.-2.*mem)*(1.+mem)/m3;

    //pipi cutoff
    double pi2phase = PKinematics::pcms(mass,makeStaticData()->GetParticleMass("pi0"),
					makeStaticData()->GetParticleMass("pi0"));

    if ((wid>0) && (use_pi_cutoff)) {
	*width= wid*pow(pi2phase,3./2.) ;
    } else if ((wid>0) && (use_hadronic_ps)) {
	Double_t pole_width = 
	    makeDynamicData()->GetParticleTotalWidthSum(
		makeStaticData()->GetParticleMass(parent_id),parent_id,1);
	Double_t mass_width = 
	    makeDynamicData()->GetParticleTotalWidthSum(
		mass,parent_id,1);

	*width= wid*(mass_width/pole_width);
    } else if (wid>0) {
	*width= wid;
    } else *width=0;


    return kTRUE;
}
#endif

Double_t PEEDirectDecay::GetWeight(Double_t *mass, Int_t *didx) {
    Warning("GetWeight","not implemented");
    return 0;
}
 

Double_t PEEDirectDecay::Eval(Double_t x, Double_t y , Double_t z , Double_t t ) const
{
    Double_t res;
    Double_t mass=x;

    if (draw_option==0) {
	return ((PChannelModel*)this)->GetWeight(&mass);
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
