//TITLE <b>Quasi-free reactions:</b> A (dummy) model for fermi samping inside a nucleus

//This class is here only for didactical reasons. It has no real physics content
//See nucleus_fermi plugin for a real implementation

#include "../src/PChannelModel.h"

//Class definition

class PFermiNucleus : public PChannelModel  {
  
 public:

    using PDistribution::GetWeight;
    PFermiNucleus(Char_t *id, Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;
    Bool_t Init(void);
    Bool_t SampleMass(void);
    Bool_t SampleMomentum(void);

 private:

    double SampleFermi(double & px, double & py, double & pz);
    PParticle * beam;
    PParticle * target;
    PParticle * spectator;
    PParticle * parent;
    PParticle * p2,*participant;
    PParticle * composite;
    
    ClassDef(PFermiNucleus,0)  //
};

PDistribution* PFermiNucleus::Clone(const char*) const {
    //clone the object
    return new PFermiNucleus((const PFermiNucleus &)* this);
};

PFermiNucleus::PFermiNucleus(Char_t *id, Char_t *de, Int_t key) : PChannelModel(id, de,key) {
    //Constructor
    beam = NULL;
    target = NULL;
    spectator = NULL;
    parent = NULL;
    participant = NULL;
    p2 = NULL;
    composite = NULL;
} ;


double PFermiNucleus:: SampleFermi(double & px, double & py, double & pz) {
    // dummy function
  
    double p = 0.5; //150 MeV
    //double p = 0.;
    double theta = acos(1.-2.*PUtils::sampleFlat());
    double phi = 2.*TMath::Pi()*PUtils::sampleFlat();
    double sth=sin(theta);
    px = p*cos(phi)*sth;
    py = p*sin(phi)*sth;
    pz = p*cos(theta);
    return p;
}

Bool_t PFermiNucleus:: Init(void) {
    
    beam = GetParticle("beam");
    target = GetParticle("target");
    parent=GetParticle("parent");

    if (!beam || !target) {
	Warning("Init","beam or target not found");
	return kFALSE;
    }

    spectator= GetParticle("spectator");
    participant=GetParticle("participant");
    p2=GetParticle("p2");

    composite= GetParticle("composite");

    return kTRUE;
}

Bool_t PFermiNucleus:: SampleMass(void) {

    Double_t massS, eS, eP, ptot, px, py, pz, t=-1., mtarget;
    while (t<0.) {
	ptot = SampleFermi(px,py,pz);                  // Fermi momentum
	massS = spectator->M();                        // mass of spectator nucleon
	eS = sqrt(ptot*ptot + massS*massS);            // spectator total energy in deuteron c.m.
	mtarget = target->M();
	t = pow(mtarget-massS,2) - 2.*mtarget*(eS-massS);  // off-shell mass**2 of participant
    }

    eP = sqrt(ptot*ptot + t);         // participant total energy

    participant->SetPxPyPzE(-px,-py,-pz,eP);  // initialize participant nucleon
    spectator->SetPxPyPzE(px,py,pz,eS);      // initialize spectator nucleon

//     participant->Print();
//     spectator->Print();

    participant->Boost(target->BoostVector()); 
    spectator->Boost(target->BoostVector());  

    *p2=*beam;

    //go into parent frame
    participant->Boost(-parent->BoostVector());
    p2->Boost(-parent->BoostVector());
    spectator->Boost(-parent->BoostVector());

    composite->Reconstruct(); //reset mass after p1 and p2 have been setted

    //boost scatter back to lab
    participant->Boost(parent->BoostVector());
    p2->Boost(parent->BoostVector());

    spectator->SetW(parent->W());                  // copy parent weight to spectator
    
    return kTRUE;

}

Bool_t PFermiNucleus:: SampleMomentum(void) {
    return kTRUE;
}


ClassImp(PFermiNucleus)


