/**************************************************
**      This class is licensed under LGPL        **
**  - free for personal and commercional use -   **
**                                               **
**   ---------------------------------------     **
**  Quasi-free scattering in the g+A scattering  **
**   ---------------------------------------     **
**                                               **
**   Author: L. Witthauer & M. Dieterle 2009     **
**                  Version 1.0                  **
**                                               **
**************************************************/

#include "PFermiMomentumGA.h"
#include "PDynamicData.h"

ClassImp(PFermiMomentumGA)
    PFermiMomentumGA::PFermiMomentumGA() {

} ;


PFermiMomentumGA::PFermiMomentumGA(const Char_t *id, const Char_t *de, Int_t key) : 
    PChannelModel(id, de, key) {
    
    fermi_model = NULL;    
    beam = NULL;
    target = NULL;
    spectator = NULL;
    parent = NULL;
    participant = NULL;
    p2 = NULL;
    composite = NULL;
    relative_warning=0;
}

PDistribution* PFermiMomentumGA::Clone(const char*delme) const {
return new PFermiMomentumGA((const PFermiMomentumGA &)* this);
};


//-----------------------------------------------------------------------------------

Bool_t PFermiMomentumGA::Init(void) {

    beam   = GetParticle("beam");
    target = GetParticle("target");
    parent = GetParticle("parent");

    if (!beam || !target) {
	Warning("Init","<%s> beam or target not found",GetIdentifier());
	return kFALSE;
    }

    spectator   = GetParticle("spectator");
    participant = GetParticle("participant");
             p2 = GetParticle("p2");

    composite= GetParticle("composite");



    return kTRUE;
}


//------- Choose arbitrary fermi momentum according to distributions in PFermiDistribution ------

Double_t PFermiMomentumGA::GetRandomFermiMomentum( Double_t & px, Double_t & py, Double_t & pz){

    //BUGBUG: What about the cross section
    if (!fermi_model) {
	//Read secondary fermi model
	//This can only happen in Init, so after the freeze-out of the user
	Int_t parent_id=makeStaticData()->GetDecayParentByKey(primary_key);
	
	Int_t id1 = parent_id/1000;
	Int_t id2 = parent_id%1000;
	if (id1 > 600) {
	    fermi_model = makeDynamicData()->
		GetParticleSecondaryModel(makeStaticData()->GetParticleName(id1), "fermi");
	} else if (id2 > 600) {
	    fermi_model = makeDynamicData()->
		GetParticleSecondaryModel(makeStaticData()->GetParticleName(id2), "fermi");
	} else if (makeStaticData()->GetParticleBaryon(id1)
		   > makeStaticData()->GetParticleBaryon(id2))  {
	    fermi_model = makeDynamicData()->
		GetParticleSecondaryModel(makeStaticData()->GetParticleName(id1), "fermi");
	} else if (makeStaticData()->GetParticleBaryon(id2)
		   > makeStaticData()->GetParticleBaryon(id1))  {
	    fermi_model = makeDynamicData()->
		GetParticleSecondaryModel(makeStaticData()->GetParticleName(id2), "fermi");
	} 

	if (!fermi_model) {
	    Error("GetRandomFermiMomentum","No fermi model found");
	    return 0.;
	}
    }

    Double_t p = (fermi_model->GetRandom());
    Double_t theta = acos(1.-2.*PUtils::sampleFlat());
    Double_t phi   = 2.*TMath::Pi()*PUtils::sampleFlat();
    Double_t sth   =sin(theta);

    px = p*cos(phi)*sth;
    py = p*sin(phi)*sth;
    pz = p*cos(theta);

    return p;

}

Bool_t PFermiMomentumGA::SampleMass(void) {
    
    Double_t massS, eS, eP, ptot, px, py, pz, t=-1., mAi;

    parent->Reconstruct(); 
 
    //    cout << "PFermiMomentumGA::SampleMass" << endl;
    while (t<0.) {
	ptot  = GetRandomFermiMomentum(px,py,pz);         // Fermi momentum
	massS = spectator->M();                           // mass of spectator fragment
	eS    = sqrt(ptot*ptot + massS*massS);            // spectator total energy in nucleus c.m.
	mAi   = target->M();
	t     = pow(mAi-massS,2) - 2.*mAi*(eS-massS);     // off-shell mass**2 of participant
    }

    eP = sqrt(ptot*ptot + t);                             // participant total energy

    participant->SetPxPyPzE(-px,-py,-pz,eP);              // initialize participant 
    spectator->SetPxPyPzE(px,py,pz,eS);                   // initialize spectator fragment
    
    //participant->Print();

    participant->Boost(target->BoostVector());
    spectator->Boost(target->BoostVector());  

    *p2=*beam;

    //    beam->Print();

    //go into parent frame
    participant->Boost(-parent->BoostVector());
    p2->Boost(-parent->BoostVector());
    spectator->Boost(-parent->BoostVector());

    composite->Reconstruct();                        

   // cout << "*********" << endl;
   //  composite->Print();
   //  spectator->Print();

    //boost scatter back to lab
    participant->Boost(parent->BoostVector());
    p2->Boost(parent->BoostVector());

    spectator->SetW(parent->W());                         // copy parent weight to spectator fragment

    return kTRUE;
}

Bool_t PFermiMomentumGA::SampleMomentum(void) {
    return kTRUE;
}


