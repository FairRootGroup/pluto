//TITLE <UserClass> Very small acceptance filter (setting variables for batch) including opening angle.


#include "../src/PBulkInterface.h"
#include "../src/PUtils.h"

class PHadesAcc : public PBulkInterface {

private:
    Double_t *hadacc, *opang;

 public:
    
    PHadesAcc();

    Bool_t Modify(PParticle ** stack, int *decay_done, int * num, int stacksize);  //bulk interface

    ClassDef(PHadesAcc,0) 
};

PHadesAcc::PHadesAcc() {
    hadacc = makeStaticData()->GetBatchValue("_hadacc");
    opang  = makeStaticData()->GetBatchValue("_opang");
}

bool PHadesAcc::Modify(PParticle ** stack, int *decay_done, int * num, int stacksize) {

    *hadacc=1.;

    for (int i=0; i< *num; i++) {
	PParticle * cur = stack[i];
//	cur->Print();




	if (cur->Is("pi+") || cur->Is("pi-")|| cur->Is("p") 
	    || cur->LeptonN()) {

#if 0
	    double momentum_resolution = 0.05;
	    double angular_resolution = 0.1*TMath::Pi()/180.;
//	    double angular_resolution = 0;

	    double mom = cur->P();

	    double theta = cur->Theta();
	    double phi   = cur->Phi();

	    double mom_measured   = PUtils::sampleGaus(mom, mom*momentum_resolution);
	    double theta_measured = PUtils::sampleGaus(theta, angular_resolution);
	    double phi_measured   = PUtils::sampleGaus(phi, angular_resolution);

	    cur->SetRho(mom_measured);
	    cur->SetTheta(theta_measured);
	    cur->SetPhi(phi_measured);

 	    cur->ResetE();
#endif

	    

	    if (cur->Theta() < ((18./180)*TMath::Pi())){
		*hadacc=0.;
		//cur->Print();
	    }
	    if (cur->Theta() > ((88./180)*TMath::Pi())){
		*hadacc=0.;
		
	    }
	    if (cur->P() < 0.1)  {
		*hadacc=0.;
		
	    }

	}

    }

    *opang=1.;

    //the opang:
    for (int i=0; i< *num; i++) {
	PParticle * ep = stack[i];
	if (ep->Is("e+")) {
	    for (int i=0; i< *num; i++) {
		PParticle * em = stack[i];
		if (em->Is("e-")) {
		    Double_t oa= ep->Vect().Angle(em->Vect());
		    if (oa < ((4./180)*TMath::Pi())){
			*opang=0.;
		    }
		}
	    }
	    
	}
    }


    return kTRUE;
};

















