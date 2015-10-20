/////////////////////////////////////////////////////////////////////
//
// Quasi-free scattering in the d+d scattering
//
//                                  Author: I. Froehlich
//
/////////////////////////////////////////////////////////////////////

#include "PFermiMomentumDD.h"

PFermiMomentumDD::PFermiMomentumDD(Char_t *id, Char_t *de, Int_t key) : 
    PFermiMomentum(id, de, key) {

}

PDistribution* PFermiMomentumDD::Clone(const char*delme) const {
    return new PFermiMomentumDD((const PFermiMomentumDD &)* this);
};

Bool_t PFermiMomentumDD:: Init(void) {

    beam = GetParticle("beam");
    target = GetParticle("target");
    parent=GetParticle("parent");

    if (!beam || !target) {
	Warning("Init","<%s> beam or target not found",GetIdentifier());
	return kFALSE;
    }

    spectator1= GetParticle("spectator");
    spectator2= GetParticle("spectator");
    p1=GetParticle("p1");
    p2=GetParticle("p2");

    composite= GetParticle("composite");
    if (!composite) {
	Warning("Init","<%s> composite not found",GetIdentifier());
	return kFALSE;
    }

    return kTRUE;
}


Bool_t PFermiMomentumDD:: SampleMass(void) {
   
    PParticle participant1, participant2;
    Int_t pair1_is_beam = 0;

    //First we have to find out which kind of reaction we have
    //for pp and nn we have no ambiguity
    if ((spectator1->ID() == 14) && (spectator2->ID() == 14)){
	participant1.SetID(13);
	participant2.SetID(13);
	pair1_is_beam = 1; //dont care
    } else if ((spectator1->ID() == 13) && (spectator2->ID() == 13)){
	participant1.SetID(14);
	participant2.SetID(14);
	pair1_is_beam = 1; //dont care
    } else if ((spectator1->ID() == 14) && (spectator2->ID() == 13)){
	//spectator1 is proton.
	//We have to check if "beam" contains the scattered proton:
	Double_t dummy = 0;
	if (beam->GetValue(P_SCATTER, &dummy)) //value used
	    pair1_is_beam = 1; //p from beam
	else if (!target->GetValue(P_SCATTER, &dummy)) {
	    //randomize if value not used
	    pair1_is_beam =(PUtils::sampleFlat() < 0.5 ? 0 : 1);
	}	
	participant1.SetID(13);
	participant2.SetID(14);
    } else if ((spectator1->ID() == 13) && (spectator2->ID() == 14)){
	//spectator2 is proton.
	//We have to check if "target" contains the scattered proton:
	Double_t dummy = 0;
	pair1_is_beam = 1; 
	if (target->GetValue(P_SCATTER, &dummy)) //value used
	    pair1_is_beam = 0; //p from target
	else if (!beam->GetValue(P_SCATTER, &dummy)) {
	    //randomize if value not used
	    pair1_is_beam =(PUtils::sampleFlat() < 0.5 ? 0 : 1);
	}

	participant1.SetID(14);
	participant2.SetID(13);
    } else {
	Warning("SampleMass","Unknown reaction");
    }
    
    Double_t massS, eS, eP, ptot, px, py, pz, t=-1., mdeut;

    //sample participant1
    while (t<0.) {
	ptot = SampleFermi(px,py,pz);                  // Fermi momentum
	massS = spectator1->M();                       // mass of spectator nucleon
	eS = sqrt(ptot*ptot + massS*massS);            // spectator total energy in deuteron c.m.
	mdeut = makeStaticData()->GetParticleMass("d");
	t = pow(mdeut-massS,2) - 2.*mdeut*(eS-massS);  // off-shell mass**2 of participant
    }

    if (t<0.) {
	Warning("SampleMass","Insufficient energy");
	return kTRUE;
    }
    eP = sqrt(ptot*ptot + t);         // participant total energy

    participant1.SetPxPyPzE(-px,-py,-pz,eP);  // initialize participant nucleon
    spectator1->SetPxPyPzE(px,py,pz,eS);      // initialize spectator nucleon

     //sample participant2
    while (t<0.) {
	ptot = SampleFermi(px,py,pz);                  // Fermi momentum
	massS = spectator2->M();                       // mass of spectator nucleon
	eS = sqrt(ptot*ptot + massS*massS);            // spectator total energy in deuteron c.m.
	mdeut = makeStaticData()->GetParticleMass("d");
	t = pow(mdeut-massS,2) - 2.*mdeut*(eS-massS);  // off-shell mass**2 of participant
    }

    if (t<0.) {
	Warning("SampleMass","Insufficient energy");
	return kTRUE;
    }
    eP = sqrt(ptot*ptot + t);         // participant total energy

    participant2.SetPxPyPzE(-px,-py,-pz,eP);  // initialize participant nucleon
    spectator2->SetPxPyPzE(px,py,pz,eS);      // initialize spectator nucleon   

    //Up to now we are in the (breakup-)deuteron frame.
    //Let us go into the lab frame first

    if (pair1_is_beam) {
	//pair 1 is beam
	participant1.Boost(beam->BoostVector()); //go from Lab to D frame
	spectator1->Boost(beam->BoostVector());  //go from Lab to D frame
	participant2.Boost(target->BoostVector()); //go from Lab to D frame
	spectator2->Boost(target->BoostVector());  //go from Lab to D frame
	if (participant2.ID() == p1->ID()) { //Identify the scattered nucleon
	    *p1=participant2;
	    *p2=participant1;
	}
	else  { //Identify the scattered nucleon
	    *p2=participant1;
	    *p1=participant2;
	}
    } else {
	//pair 2 is beam
	participant2.Boost(beam->BoostVector()); //go from Lab to D frame
	spectator2->Boost(beam->BoostVector());  //go from Lab to D frame
	participant1.Boost(target->BoostVector()); //go from Lab to D frame
	spectator1->Boost(target->BoostVector());  //go from Lab to D frame
	if (participant2.ID() == p1->ID()) { //Identify the scattered nucleon
	    *p1=participant2;
	    *p2=participant1;
	}
	else  { //Identify the scattered nucleon
	    *p2=participant1;
	    *p1=participant2;
	}
    }


    //go into parent frame
    p1->Boost(-parent->BoostVector());
    p2->Boost(-parent->BoostVector());
    spectator1->Boost(-parent->BoostVector());
    spectator2->Boost(-parent->BoostVector());

    composite->Reconstruct(); //reset mass after p1 and p2 have been setted

    //boost scatter back to lab
    p1->Boost(parent->BoostVector());
    p2->Boost(parent->BoostVector());

    spectator1->SetW(parent->W());                  // copy parent weight to spectator
    spectator2->SetW(parent->W());

    return kTRUE;
}


ClassImp(PFermiMomentumDD)
