////////////////////////////////////////////////////////
//  Pluto interface to FairRoots primary generator
//
//
//                    Author:  Ingo Froehlich
//                    Written: 13/10/2009
//
////////////////////////////////////////////////////////

#include "PFairGenerator.h"
#include "PDynamicData.h"
#include "PChannel.h"
#include "PDistributionManager.h"

#include "TDatabasePDG.h"

PFairGenerator::PFairGenerator() {

    fPriority=99999; //a really high number!

}

bool PFairGenerator::Modify(PParticle ** stack, int *decay_done, int * num, int stacksize) {
	    
    fNumberParticles = 0;

    //cout << "PFairGenerator::Modify" << endl;

    for (int i=0;i<*num;i++) {

	//Loop over stack particles
	if (stack[i]->IsActive()) { 
	    //count only active particles
	    //because the fireball macros have a multiplicity sampling
	    if (!decay_done[i]) {
		//count only undecayed particles

		if (fNumberParticles == FAIRGENERATOR_STACKSIZE) {
		    Error("Modify","FAIRGENERATOR_STACKSIZE reached");
		    return kFALSE;
		}
		fLocalStack[fNumberParticles]=stack[i];
		fNumberParticles++;
		
	    }
	}
    }

    return kTRUE;
}

Bool_t PFairGenerator::GetNextParticle(Int_t *pdgType, Double_t *px,  Double_t *py,  
				       Double_t *pz,  Double_t *vx,  Double_t *vy,  Double_t *vz) {
    //iterator for reading the particles

    if (!fNumberParticles) return kFALSE;

    TDatabasePDG* dataBase = TDatabasePDG::Instance();

    fNumberParticles--;
    PParticle *part = fLocalStack[fNumberParticles];
    
    TLorentzVector mom = part->Vect4();
    *px = mom.Px();
    *py = mom.Py();
    *pz = mom.Pz();

    TVector3 vertex = part->getVertex();
    *vx = vertex.x();
    *vy = vertex.y();
    *vz = vertex.z();

    *pdgType = dataBase->ConvertGeant3ToPdg( part->ID() );

 

    return kTRUE;

}


ClassImp(PFairGenerator) 
