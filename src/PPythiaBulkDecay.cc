////////////////////////////////////////////////////////
//  Pythiao bulk decay base 
//
//  This class let all particles with tau<tauMax 
//  decay in the Pythia way
//
//  BUGBUG: Class is UNTESTED
//
//                    Author:  Ingo Froehlich
//                    Written: 10/07/2007
//                    Revised:
//
////////////////////////////////////////////////////////

#include "PPythiaBulkDecay.h"
#include "PDynamicData.h"
#include "PChannel.h"
#include "PDistributionManager.h"

#ifdef USE_PYTHIA6
  #include "TPythia6.h"
#endif


PPythiaBulkDecay::PPythiaBulkDecay() {
#ifndef USE_PYTHIA6
    Warning("PPythiaBulkDecay","Pythia library not included");
#endif

    fPythia=NULL;

}

bool PPythiaBulkDecay::Modify(PParticle ** stack, int *decay_done, int * num, int stacksize) {
    //decay all particle
    //input: particle array with *num members
    //new particles should be instantiated and *num increased

#ifndef USE_PYTHIA6
    Warning("PPythiaBulkDecay","Pythia library not included");
    return kFALSE;
#endif 
    
    int maxnp=stacksize-*num;
    if (!fPythia) return kFALSE;
    
//    Int_t child_ids[maxnp+1];
    PParticle *cur_p, *work[maxnp+1];
    int st_i1 = 0;
#ifdef USE_PYTHIA6
    int st_i2 = *num;
#endif 
    int st_i3 = *num;
    while (st_i1 < st_i3) {
        cur_p = work[0] = stack[st_i1];
	
	
	if (cur_p->IsActive() && 
	    makeDynamicData()->GetParticleLife(cur_p->ID()) 
	    < tauMax && !decay_done[st_i1]) {

#ifdef USE_PYTHIA6
            fPythia->SetN(1);      // load particle into Pythia
            fPythia->SetK(1,1,1);  // status code  -> manual p. 54 
            fPythia->SetK(1,2,PData::idtokf(cur_p->ID()));  // id 
            fPythia->SetK(1,3,0);
            fPythia->SetK(1,4,0);
            fPythia->SetK(1,5,0);
            fPythia->SetP(1,1,cur_p->Px()); // px in GeV/c
            fPythia->SetP(1,2,cur_p->Py()); // py
            fPythia->SetP(1,3,cur_p->Pz()); // pz
            fPythia->SetP(1,4,cur_p->E());  // E
            fPythia->SetP(1,5,cur_p->M());  // mass
            fPythia->SetV(1,1,cur_p->X());  // creation vertex in mm
            fPythia->SetV(1,2,cur_p->Y());
            fPythia->SetV(1,3,cur_p->Z());
            fPythia->SetV(1,4,cur_p->T());  // production time in mm/c
            fPythia->SetV(1,5,cur_p->getProperTime()); // proper time in mm/c

            fPythia->Pyexec();  // decay particle

            Int_t np = fPythia->GetN()-1;   // number of decay products
            if (np > maxnp) {
		Warning("Modify","np=%d > maxnp=%d\n",np,maxnp);
		np = maxnp;
            }
            if (np>1) cur_p->setDaughterIndex(st_i3);
	    
            Int_t np1st=0;   // number of 1st generation decay products
	    
            for (k=1; k<=np; k++) {  // 1st loop needed to setIndex() of all siblings
		Int_t k1=k+1;
		if (fPythia->GetK(k1,3) != 1) break;   // stop if parent is not
		// 1st in list
		Int_t prodId = PData::kftoid(fPythia->GetK(k1,2));
		if (!prodId) printf("PReaction::loop(): unknown decay product\n");
		work[k] = new(p_array[st_i3-st_i2]) PParticle(prodId);
		stack[st_i3] = work[k];
		stack[st_i3]->setIndex(st_i3+1);
		stack[st_i3]->setParent(work[0]);
		st_i3++;
		np1st++;
            }

	    for (k=1; k<=np1st; k++) {
		Int_t k1=k+1;
		work[k]->setActive(); // activate children
		work[k]->setSourceId(cur_p->getSourceId());
		work[k]->setParentId(cur_p->ID());
		work[k]->setParentIndex(cur_p->getIndex());
		work[k]->setDaughterIndex(-1);
		if (k!=np1st) {
		    work[k]->setSiblingIndex(work[k+1]->getIndex());
		    work[k]->setSibling(work[k+1]);
		}
		else {
		    work[k]->setSiblingIndex(work[1]->getIndex());
		    work[k]->setSibling(work[1]);
		}
		work[k]->setVertex(fPythia->GetV(k1,1),fPythia->GetV(k1,2),
				   fPythia->GetV(k1,3),fPythia->GetV(k1,4)); // vertex
		work[k]->setProperTime(fPythia->GetV(k1,5));                 // proper time
		work[k]->SetPxPyPzE(fPythia->GetP(k1,1),fPythia->GetP(k1,2),
				    fPythia->GetP(k1,3),fPythia->GetP(k1,4));
		work[k]->SetW(cur_p->W());
//            if (!allPARTICLES && (fPythia->GetP(k1,4)!=0)) // has daughter(s)
//                             work[k]->setInActive();
            }
	    decay_done[st_i1]=1;	    
#endif
	}

    }
    
    return kTRUE;
}

ClassImp(PPythiaBulkDecay) 
