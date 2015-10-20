/////////////////////////////////////////////////////////////////////
//  PChannel Class implementation file
//
//  PChannel is one of the main classes of Pluto.
//  A Channel object represents one step of a reaction process,
//  namely of a parent particle decaying into n decay products.
//  A Channel consists of parent and decay particles and
//  the decay mode. For an overview of the Pluto package
//  see also PData, PParticle, PReaction, and PFilter.
//
//                                  Author:  M.A. Kagarlis
//                                  Written: 21.01.99
//                         p+d part revised: 13.10.00 by V. Hejny
//                                  Revised: 19.08.03 by R. Holzmann
//            angular distributions Revised: 08.11.03 by I. Froehlich
//           propagate parent index Revised: 22.03.05 by R. Holzmann
//                     set sourceId Revised: 13.04.05 by R. Holzmann
//                  helicity angles Revised: 28.11.05 by I. Froehlich
//          eta->3pi matrix element Revised: 05.12.05 by I. Froehlich
//             ds_dt(cos) bug fixed Revised: 11.07.06 by R. Holzmann
//            e_cm static bug fixed Revised: 29.08.06 by R. Holzmann
//
//
////////////////////////////////////////////////////////////////////

using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "Pdefines.h"

const char *Message[14] = {
    "PChannel: operation successful",
    "PChannel: invalid pids",
    "PChannel: reconstruct failed, increase beam energy",
    "PChannel: decay failed, system invariant mass is zero",
    "PChannel: decay failed, insufficient energy",
    "PChannel: genbod failed, fewer than two decay products",
    "PChannel: genbod failed",
    "PChannel: baryon_cos algorithm failed",
    "PChannel: mflag must be 0 or 1",
    "PChannel: aflag must be 0 or 1",
    "PChannel: bflag must be 0 or 1",
    "PChannel: reconstruct failed for deuteron Fermi sampling",
    "PChannel: invalid decay-mode index in the constructor",
    "PChannel: inconsistent particle array with decay-mode index"};

#include "PChannel.h"
#include "PDynamicData.h"
#include "PUtils.h"
#include "PFileInput.h"
#include "PFireball.h"
#include "PDiLepton.h"
#include "PKinematics.h"



double PChannel::globalWeight = 1.;              // global channel weight

PChannel:: PChannel(PParticle **particles, int nt, int mf, int af, int bf) {
    // Channel constructor by pointer to array of pointers to particles
    // (parent, decay products), number of decay particles (default 2).
    // ____________________________________________________________________________
    // Flags: (obsolete)

    makeStdData()->fillDataBase();//init physics, if not yet done
    init_done=kFALSE;
    e_cm=0;
    int i, k=-1;
    status=(!particles);         // true if argument pid array is invalid
    n=nt;                        // number of product particles
  
    if (!status) {
	ipid.Set(n+1);                 // pid array
	for (i=0;i<=n;++i) {
	    k = particles[i]->ID();      // id of current particle
//      if ((i&&k>1000)||!k) {       // error in pid (Not longer true, IF)
	    if (!k) { 
		status=1;
		ipid.~TArrayI();
		break;
	    }
	    ipid[i]=k;
	}
    }
    if (status) {                  // failed
	printf("%s:",Message[status]);
	for(i=0;i<=n;i++) printf(" %d",particles[i]->ID());
	printf("\n");
	return;
    }

    ptcls = particles;
    ptcls_extern = -1;

    IsReaction();     // if reaction channel, identify beam, target & decay mode



    particles[0]->SetDecayModeIndex(makeStaticData()->GetDecayIdxByKey(decay_key));
    // cout << "Parent didx set to " << makeStaticData()->GetDecayIdxByKey(decay_key) << endl;
    w=particles[0]->W();           // weight on instantiation
    ecm=particles[0]->M();         // system invariant mass on instantiation
    DMIndex=makeStaticData()->GetDecayIdx(ipid.GetArray(),n); // match decay mode with PData

}

PChannel:: PChannel(int idx, PParticle ** particles, int mf, int af, int bf) {
    // Channel constructor by decay-mode index
    
    makeStdData()->fillDataBase();//init physics, if not yet done
    init_done=kFALSE;
    e_cm=0;
//    status=(idx<0||idx>PData::maxnummodes);
    //get key
    Int_t key=makeDataBase()->GetEntryInt("didx", idx);
    Int_t ia[11];
    ia[0]=10;
    if (key<0) {
	printf("%s\n",Message[12]);    // invalid decay mode
	return;
    } else makeStaticData()->GetDecayMode(idx,ia); // retrieve decay-mode info
    DMIndex=idx;                     // decay-mode index

    n=ia[0];                         // number of products in current mode

    ipid.Set(n+1);                   // set up pid array
    
    if (!particles) {                // no particle-array passed as argument
	ptcls = new PParticle*[n+1];   // create channel particle array
	for(int i=0; i<=n; i++) ptcls[i] = 0;
	ptcls_extern = 0;
    }
    else {
	ptcls = particles;          // copy
	ptcls_extern = -1;
    }
    
    int i;
    for (i=0;i<=n;++i) {             // fill particle and id arrays
	if (particles) ipid[i]=particles[i]->ID();
	else {                         // create from decay-mode info
	    if (!i)                      // parent id
		ipid[0]=makeStaticData()->GetDecayParent(idx);
	    else ipid[i]=ia[i+1];        // product ids
	    ptcls[i]=new PParticle(ipid[i]);
	}
    }      
    
    if (particles) {                 // consistency checks
	ia[0]=ipid[0];                 // reset first entry to parent pid
	PUtils::isort(ia,n+1);         // sort pids (parent & products) in ascending order
	TArrayI test(ipid);            // make a copy of channel pid array
	PUtils::isort(test.GetArray(),n+1);// sort as well
	int hit=0;
	for (i=0;i<=n;++i)             // ompare id arrays
	    hit += (test[i]==ia[i]);
	if (hit!=n) {                  // input array inconsistent with decay-mode index
	    printf("%s\n",Message[14]);
	    ipid.~TArrayI();
	    return;
	}
    }
    
    IsReaction();                  // if reaction channel identify beam, target & decay mode

    ptcls[0]->SetDecayModeIndex(makeStaticData()->GetDecayIdxByKey(decay_key));
    // cout << "Parent didx set to " << makeStaticData()->GetDecayIdxByKey(decay_key) << endl;
    w=ptcls[0]->W();               // weight on instantiation
    ecm=ptcls[0]->M();             // system invariant mass on instantiation
}     

PChannel:: PChannel(int idx, PParticle & parent, int mf, int af, int bf) {
    // Channel constructor by decay-mode index and parent reference
    
    makeStdData()->fillDataBase();//init physics, if not yet done
    init_done=kFALSE;
    e_cm=0;

    //get key
    Int_t key=makeDataBase()->GetEntryInt("didx", idx);
    Int_t ia[11];
    ia[0]=10;
    if (key<0) {
	printf("%s\n",Message[12]);    // invalid decay mode
	return;
    } else makeStaticData()->GetDecayMode(idx,ia); // retrieve decay-mode info
    DMIndex=idx;                     // decay-mode index

    n=ia[0];                         // number of products in current mode

    ipid.Set(n+1);                   // set up pid array
    
    ptcls = new PParticle*[n+1];     // create channel particle array
    ptcls[0] = &parent;
    ptcls_extern = 1;
    ipid[0] = parent.ID();           // parent id

    if (ipid[0]!=makeStaticData()->GetDecayParent(idx)) {
	printf("%s\n",Message[14]);
	ipid.~TArrayI();
	delete [] ptcls;
	return;
    }

    int i;
    for (i=1;i<=n;++i) {             // fill particle and id arrays
	ipid[i]=ia[i];                 // product ids
	ptcls[i]=new PParticle(ipid[i]);
    }

    IsReaction();                    // if reaction channel identify beam, target & decay mode

    ptcls[0]->SetDecayModeIndex(makeStaticData()->GetDecayIdxByKey(decay_key));
    // cout << "Parent didx set to " << makeStaticData()->GetDecayIdxByKey(decay_key) << endl;
    w=ptcls[0]->W();                 // weight on instantiation
    ecm=ptcls[0]->M();               // system invariant mass on instantiation
}     

void PChannel:: IsReaction() {
    // Identify beam and target if the current channel is a beam + target reaction.
 
    event_impact_param = (makeStaticData()->GetBatchValue("_event_impact_param"));
    event_plane        = (makeStaticData()->GetBatchValue("_event_plane"));
    weight_version     = makeStaticData()->GetBatchValue("_system_weight_version");

    makeDynamicData()->PushPChannel((TObject *) this);

    parent      = ptcls[0];
    orig_parent = ptcls[0];
    int k=GetParentSize();                   // get parent size
    thSrc=0;
    tcSrc=0;
    dlSrc=0;
    fEnablePattern = 0;
    thermal_disable_didx=0;
    num_not_finalized = 0;

    distribution_position = 0;
    print_tentative=kTRUE;
    quasi_pchannel=NULL;
    weight_sum=0;
    tid=bid=0;
    k=GetParentSize();
// The grandparent is not known yet at this stage, it is set in PReaction::setUp() only!

    parent->SetSourceId(parent->ID());
    parent->SetParentId(0);

    for (int i=1;i<=n;i++) {
	ptcls[i]->SetValue(CHANNEL_POS ,(double)i);
    }

    emin=0.;                       // channel energy threshold

    if (k==1) {                              // decay of size=1 parent
	if (parent->IsFireball()) thSrc=1;
	else if (parent->IsFileInput()) tcSrc=1;
	else if (parent->IsDilepton()) dlSrc=1;
	else {
	    IdChannel();

	    //get data base key
	    decay_key = makeStaticData()->GetDecayKey(ipid.GetArray(), n);
	    CheckDecayKey();
	    return;	 
	}
    }

    //Continue with compound parents

    parent->SetParentIndex(-2);              // Set some indexes for the compound object
    parent->SetIndex(-1);


    tid=ipid[0]/1000;                    // target id
    bid=ipid[0]%1000;                    // beam id

    IdChannel();                             // identify reaction if known
    //cout << "parent: " << ipid[0] << endl;
    //get data base key
    decay_key = makeStaticData()->GetDecayKey(ipid.GetArray(), n);
    if ((decay_key < 0) && (ipid[0]>999)) {
	//If not yet defined, add the decay of composites to the data base
	decay_key = makeStaticData()->AddDecay(ipid.GetArray(), n);
    }
    if (!thSrc) CheckDecayKey();
    //check for emin in initial states
    for (int i=1;i<=n;++i) { 
	emin += PData::LMass(ptcls[i]->ID());
    }
    //cout << "emin set to: " << emin << endl;
    //BUGBUG: check for an additional Decayemin

}

void PChannel::CheckDecayKey(void) const {

    if (decay_key < 0) {
	Error("Initialization","No database entry for:");
	cout << "               ";
	PrintReaction(0);
	Error("Initialization","Bye...");
	exit(1);
    }


}

void PChannel::IdChannel() {
    // Make quasi-channel for d+p reactions

    int i, nucleon_pid;

    spectator=NULL;
    participant=NULL;
    quasi_composite=NULL;
    quasi_pchannel=NULL;

    //Workaround for backward-compatibility
    //check id beam or target is a d, and we have
    //no d in the exit channel. In this case, we
    //construct a helper PChannel in oder to get the
    //"real" scattering

    //First, we check if all particles have the spec flag=-1
    //this is the case when using the decay parser but not "(x)"

    int has_spec=0;
    for (i=1;i<=n;++i) {
	if (ptcls[n]->IsSpectator() != -1)
	    has_spec=1;
    }
    
    if (!has_spec) return;

    has_spec = ptcls[n]->IsSpectator();

    //if spectator model not forced, do it only for d
    if (has_spec == 0) {
	if ((bid != makeStaticData()->GetParticleID("d")) &&
	    (tid  != makeStaticData()->GetParticleID("d")))
	    return;
	//in this case the, spectator must be a nucleon
	if (!ptcls[n]->IsNucleon()) return;
    }
        
    for (i=1;i<=n;++i) {
	if (!has_spec) {
	    if (ptcls[i]->ID() == makeStaticData()->GetParticleID("d")) 
		return; // coherent scattering
	}
	if (ptcls[i]->ID()>1000) return; //dummy PChannel
    }

    //get the nucleon pid
    if (bid == makeStaticData()->GetParticleID("d")) nucleon_pid=tid;
    else nucleon_pid=bid;

    //identify the spectator
    spectator = ptcls[n];

    //identify the participant
    if (spectator->ID() == makeStaticData()->GetParticleID("p"))
	participant=new PParticle("n");
    else if (spectator->ID() == makeStaticData()->GetParticleID("n"))
	participant=new PParticle("p");
    else {
	//loop over first part, misuse "spectator flag"
	for (i=1;i<=(n-1);++i) {
	    if (ptcls[i]->IsSpectator() == 1)
		participant=new PParticle(ptcls[i]);
	}
	//special case when we have a quasi-reaction (a b) is more complicated
	//in this case b is the participant
	Int_t a_pos=-1;
	for (i=1;i<=(n-1);++i) {
	    if (ptcls[i]->IsSpectator() == 3)
		participant=new PParticle(ptcls[i]);
	    if (ptcls[i]->IsSpectator() == 2)
		a_pos=i;
	}
	if (a_pos != 1)  {
	    Error("IdChannel","Quasi-reaction must come at first");
	}
	else if ((a_pos == 1) && (participant)) {//....  a and b have to be removed
	    for (i=1;i<=(n-3);++i) {
		ptcls[i] = ptcls[i+2];
		ipid[i] = ipid[i+2];
	    }
	    n-=2;
	}

	if (!participant) {
	    Warning("IdChannel","No valid spectator found");
	    return;
	}
    }

    PParticle dummy(nucleon_pid);
    if (nucleon_pid==tid) { //d beam
	quasi_composite = new PParticle(*participant + dummy);
    } else {
	quasi_composite = new PParticle(dummy + *participant);  
    }
    PParticle ** cc1 = new PParticle*[3];   // create channel particle array
    cc1[0]=parent;
    cc1[1]=quasi_composite;
    cc1[2]=spectator;

    quasi_pchannel = new PChannel(cc1,2);

    //now we have to touch a little bit myself:
    //the spectator is removed from the list, and the parent id replaced by
    //the dummy composite
    
    Info("IdChannel","(%s) Quasi-free production", PRINT_AUTO_ALLOC);

    n--;
  
    parent=quasi_composite;
    orig_parent=quasi_composite;
    ptcls[0]=quasi_composite;
    ipid[0]=quasi_composite->ID();
 
}

void PChannel:: ThermalSampling() {
    // parent is midrapidity thermal source

    double px, py, pz, Energy;
    float b = 0., phi = 0.;
    int i, nTherm = n;
    PFireball * pFire =          // pointer to thermal source
	(ptcls[0]->IsFireball()) ? (PFireball *)ptcls[0] : NULL;

    if (pFire->IsRandomB()) {             // sample impact parameter
	if (*event_impact_param == 0.) {
	    b = pFire->sampleB() + 0.001;     // sample new value
	    phi = 2.*TMath::Pi()*PUtils::sampleFlat(); // ev plane
	}
	else {
	    b = *event_impact_param;    // use previous value
	    phi = *event_plane;
	}
	nTherm = pFire->sampleNProd(b);
	if (nTherm > n) nTherm = n;
    }
    else if (pFire->IsRandomN()) {        // sample multiplicity directly
	nTherm = pFire->sampleNProd();
	if (*event_plane == 0.) 
	    phi = 2.*TMath::Pi()*PUtils::sampleFlat(); // ev plane
	else phi = *event_plane;
	if (nTherm > n) nTherm = n;
    }
    //    cout << "got " << nTherm << " thermal products" << endl;
    *event_impact_param = b;
    *event_plane        = phi;
    for(i=1;i<=nTherm;i++) {
	int didx = ptcls[i]->GetDecayModeIndex(1);
	//	cout << didx << endl;

	if (nTherm==1) {
	    //cout << didx << endl;
	    //didx = -1;
	    pFire->samplePartCM(px, py, pz, Energy,didx);
	    //cout << "EXIT" << endl;
	}
	else {
	    //when having more then one product with different didx
	    //it causes endless time because the 2dim ROOT sampling
	    //will regenerate the polynomials each time
	    if (thermal_disable_didx)
		pFire->samplePartCM(px, py, pz, Energy,-1);
	    else if ((i>1) && didx != ptcls[i-1]->GetDecayModeIndex(1)) {
		//inconsistent didx -> disable
		thermal_disable_didx = 1;
		Warning("ThermalSampling","Different decays, disable partial decay sampling");
	    } else {
		pFire->samplePartCM(px, py, pz, Energy,didx);
	    }
	}
	ptcls[i]->SetPxPyPzE(px, py, pz, Energy);
	ptcls[i]->SetW(w*(pFire->mtScale(ptcls[i]->M())));  // weight * mt-scaling
	ptcls[i]->SetSourceId(pFire->ID());
	ptcls[i]->SetParentId(pFire->ID());
	ptcls[i]->SetActive();
	// cout << "active " << i << endl;
    }
    for(i=nTherm+1;i<=n;i++) {
	// cout << "inactive " << i << endl;
	ptcls[i]->SetInActive();  // inactivate rest
    }

    TVector3 Beta=pFire->BoostVector();
    if (Beta.Mag()>0.) {
	//cout << "boosting" << endl;
	for(i=1;i<=nTherm;++i) {
	    //	    ptcls[i]->Print();
	    ptcls[i]->Boost(Beta);
	}
    }



}

void PChannel:: MakeDilepton() {
    // parent is dilepton source

    double px, py, pz, Energy;
    PDiLepton * pDL =          // pointer to dilepton
	(ptcls[0]->IsDilepton()) ? (PDiLepton *)ptcls[0] : NULL;

    pDL->samplePartCM(px, py, pz, Energy);
    ptcls[1]->SetPxPyPzE(px, py, pz, Energy);
    ptcls[1]->SetW(w);
    ptcls[1]->SetSourceId(pDL->ID());
    ptcls[1]->SetParentId(pDL->ID());
    ptcls[1]->SetActive();

    TVector3 Beta=pDL->BoostVector();
    if (Beta.Mag()>0.) ptcls[1]->Boost(Beta);
}


Int_t PChannel:: ReadFileInput() {  // parent is file input

    Double_t px, py, pz, Energy, vx, vy, vz, vt;
    Float_t b = 0., phi = 0.;
    Int_t id, idSrc, idPar, indPar=-2, ret;
    Int_t i, nFile = 0, nMax;
    PFileInput *pFile =             // pointer to file input interface
	(ptcls[0]->IsFileInput()) ? (PFileInput*)ptcls[0] : NULL;

    nMax = nFile = pFile->readEventHeader(b);  // read event header
    if (nFile==-1) {
	printf("**** EOF reached on file input ****\n");
	for(i=1;i<=n;i++) ptcls[i]->SetInActive();  // inactivate all
	return 8;
    }

    if (nMax > n) nMax = n;
  
    if (*event_impact_param == 0.) {
	phi = 2.*TMath::Pi()*PUtils::sampleFlat(); // ev plane
    }
    else {
	b =   *event_impact_param;           // use previous value
	phi = *event_plane;
    }

    *event_impact_param = b;
    *event_plane = phi;
    
    for(i=1;i<=nFile;i++) {
	ret=pFile->readParticle(px,py,pz,Energy,vx,vy,vz,vt,id,idSrc,idPar,indPar,w);
	if (ret==-1) break;  // EOF
	if (i>n) continue;
	ptcls[i]->SetID(id);           // reset particle ID
	ptcls[i]->SetPxPyPzE(px, py, pz, Energy);
	ptcls[i]->ResetE();            // make Energy consistent with mass
	ptcls[i]->SetW(w);
	ptcls[i]->SetSourceId(idSrc);
	ptcls[i]->SetParentId(idPar);
	ptcls[i]->SetParentIndex(indPar);
	ptcls[i]->SetActive();
	ptcls[i]->SetVertex(vx,vy,vz,vt);
    }
    if (ret==-1) return 8;
    for(i=nMax+1;i<=n;i++) ptcls[i]->SetInActive();  // inactivate rest

    TVector3 Beta=pFile->BoostVector();
    if (Beta.Mag()>0.) for(i=1;i<=nMax;++i) ptcls[i]->Boost(Beta);
    return 0;
}

int PChannel:: CheckSiblings(PParticle * p, PDistribution * dist, int flag) {
    if (p == NULL) {
	return -1;
    }

    if (p->GetParent()) //Fireball daughters: no check for siblings
	if ((p->GetParent()->ID() > 500) && (p->GetParent()->ID() <1000)) return -1;

    if (p->ID() >= 1000) return -1; //no sibling for compounds

   //  cout << "******************checkSiblings for " <<dist->GetName() <<   ", SEED: " << endl;
//      p->Print();
    
    PParticle * currentSibling = p->GetSibling();   
    Int_t counter=0; //break condition
    Int_t ret    =-1;
    while (currentSibling != p) {
	if (currentSibling == NULL) return ret;
	if (ret == 0) { //found one more "additional" sibling -> reset 
	    ret = -1;
	    //cout << "RESET" << endl;
	}
  	//cout << "My sibling is:" << endl;
 	//currentSibling->Print();

	if (dist->SetParticle(currentSibling, currentSibling->ID(), flag | PARTICLE_LIST_SIBLING) == 0) { //I found the missing sibling
	    ret = 0;
	    //cout << "found the missing sibling" << endl;
	}  
	//jump to next sibling
	currentSibling = currentSibling->GetSibling();
	counter++;
	if (counter==8) { //max 7 decay particles
	    Warning("CheckSiblings","Sibling chain seems to be corrupted");
	    Print();
	    return -1;
	}
    }
//     cout << "RET: " << ret << endl;
//     dist->Print();
    if (ret == -1) dist->ResetRelatives(flag | PARTICLE_LIST_SIBLING); //RESET current flag
    return ret;
}

Bool_t PChannel:: Reset() {
    for (int j=0;j<distribution_position;j++) {
	dist[j]->ResetRelatives();
	dist[j]->SetEnable(1); //can be done as disabled distr. are never attached
    }

    grandparent=NULL;
    grandgrandparent=NULL;
    init_done=kFALSE;

    return kTRUE;
}

Bool_t PChannel::SetDaughters() {
    if (thSrc)  return kTRUE; //makes no sense and print errors
    parent->ResetDaughters();
    
    if (parent->IsActive() == kTRUE) {
	for (Int_t k=1;k<=n;k++) {
	    parent->SetDaughter(ptcls[k]);
	}
    }
    return kTRUE;
}

Bool_t PChannel::Init() {
    //On the 1st call, set the grandparent, etc...
    //look for tentative distributions
    
    int has_printed = 0;
    //set parent, grandparent
    if (orig_parent!=ptcls[0]) {
	//Workaround for PlutoBulkDecay
	parent = ptcls[0];
	*orig_parent=*(ptcls[0]);
    }

    grandparent=parent->GetParent();
    grandgrandparent=NULL;
    if (grandparent) grandgrandparent=grandparent->GetParent();

   

    Double_t set_mode_index = 0;

    for (int j=0;j<distribution_position;j++) { //loop over all valid distributions
	if (dist[j]->GetVersionFlag() & VERSION_INVERT_WEIGHTING) {
	    set_mode_index++;
	} else if ((dist[j]->GetVersionFlag() & VERSION_WEIGHTING) 
		   && (dist[j]->GetVersionFlag() & VERSION_GENERATOR_MC) 
//		   && (dist[j]->GetExpectedWeightMean() < 0) 

// This must be disabled, otherwise the Delta mass is not taken correctly
// Mass shape of the epja figure disagrees for sampling/weighting method		   
	    ) {
	    set_mode_index++;
	}
	
    }

    if (set_mode_index < 2)  //
	ptcls[0]->SetDecayModeIndex(makeStaticData()->GetDecayIdxByKey(decay_key));
    else {
// 	cout << "partial sampling disabled" << endl;
// 	Print();
	ptcls[0]->SetDecayModeIndex(-1);
    }


    if (parent->IsActive() == kTRUE) {
	parent->SetDaughterIndex(ptcls[1]->GetIndex());
	Int_t parId, parInd;
	for (Int_t k=1;k<=n;k++) {
	    ptcls[k]->SetActive();          // activate children
	    ptcls[k]->SetSourceId(parent->GetSourceId());


	    parId = parent->ID();
	    parInd = parent->GetIndex();
	    if ((parId==50 || parId==51) && grandparent!=NULL) {
		parId += 1000*(grandparent->ID()); // if parent = dilepton,
		// add grandparent id
		parInd = grandparent->GetIndex();
	    }
	    ptcls[k]->SetParentId(parId);
	    ptcls[k]->SetParentIndex(parInd);
	    ptcls[k]->SetDaughterIndex(-1);
	    if (k!=n) {
		ptcls[k]->SetSiblingIndex(ptcls[k+1]->GetIndex());
		ptcls[k]->SetSibling(ptcls[k+1]);
	    }
	    else {
		ptcls[k]->SetSiblingIndex(ptcls[1]->GetIndex());
		ptcls[k]->SetSibling(ptcls[1]);
	    }
	  
	}
    } //end parent active


    // check tentative distributions for additional particles

    for (int j=0;j<distribution_position;j++) {
	if ((dist[j]->GetStatus() == -1) && dist[j]->GetEnable()){ //not all particles set
	    if (has_printed == 0 && print_tentative) {
		has_printed = 1;
		cout << "checking tentative distributions for ";
		if (grandparent) {
		    if (grandparent->ID() < 1000)
			printf("%s --> ",makeStaticData()->GetParticleName(grandparent->ID()));
		}
		PrintReaction();
		cout << endl;
	    }
	    if (print_tentative) 
		printf("                             [checking] %s\n",
		       dist[j]->GetDescription());
	    //check grand(grand)parent itself

	    if (grandparent) {
		if (grandparent->ID() < 1000)
		    dist[j]->SetParticle(grandparent, grandparent->ID(), PARTICLE_LIST_GRANDPARENT);
		else
		    dist[j]->SetParticle(grandparent, DISTRIBUTION_QUASI_PID, PARTICLE_LIST_GRANDPARENT);
	    }
	    if (grandgrandparent) {
		if (grandgrandparent->ID() < 1000)
		    dist[j]->SetParticle(grandgrandparent, grandgrandparent->ID(), PARTICLE_LIST_GRANDGRANDPARENT);
		else
		    dist[j]->SetParticle(grandgrandparent, DISTRIBUTION_QUASI_PID, PARTICLE_LIST_GRANDGRANDPARENT);
	    }
#if 1
	    //sometimes we have to jump over a compound particle
	    if (parent->ID() >= 1000) {
		if (parent->GetScattering(0)) {
		    dist[j]->SetParticle(parent->GetScattering(0), parent->GetScattering(0)->ID(), PARTICLE_LIST_GRANDPARENT);
		    dist[j]->SetParticle(parent->GetScattering(1), parent->GetScattering(1)->ID(), PARTICLE_LIST_GRANDPARENT);
		}
	    } else if (grandparent) {
		if (grandparent->ID() >= 1000) {
		    if (grandparent->GetScattering(0)) {
			dist[j]->SetParticle(grandparent->GetScattering(0), grandparent->GetScattering(0)->ID(), PARTICLE_LIST_GRANDGRANDPARENT);
			dist[j]->SetParticle(grandparent->GetScattering(1), grandparent->GetScattering(1)->ID(), PARTICLE_LIST_GRANDGRANDPARENT);
		    }
		}
	    }

#endif

	    //check siblings for parent, grandparent, grandgrandparent	  
	    CheckSiblings(parent, dist[j], PARTICLE_LIST_PARENT);
	    if (grandparent) {
		CheckSiblings(grandparent, dist[j], PARTICLE_LIST_GRANDPARENT);
		if (grandgrandparent)
		    CheckSiblings(grandgrandparent, dist[j], PARTICLE_LIST_GRANDGRANDPARENT);
	    }

	    //check for granddaughters...
	    for (Int_t k=1;k<=n;k++) {
	      //looping over daughters
	      int kk=0;

	      while (ptcls[k]->GetDaughter(kk) != NULL) {  

		dist[j]->SetParticle(ptcls[k]->GetDaughter(kk), 
				     ptcls[k]->GetDaughter(kk)->ID(), 
				     PARTICLE_LIST_GRANDDAUGHTER);
		
	      
		kk++;
	      }
	      
	    }
	    


	    //final check
	    if (dist[j]->GetStatus() == -1) { //not all particles set
		if (print_tentative) printf("                             [removed] %s\n",dist[j]->GetDescription());
		dist[j]->SetEnable(0); 
	    }
	    else {
		if (dist[j]->Init() == kTRUE) {
		    if (print_tentative) printf("                             [enabled] %s\n",dist[j]->GetDescription()); }
		else {
		    printf("                             [error] %s:%s\n",dist[j]->GetDescription(),dist[j]->GetIdentifier());
		    dist[j]->SetEnable(0);
		  
		}
	    }
	}
      
    }

    //after all, look for artificial models (primary)
    //this can happen when parent is composite

    Int_t num_mass_sampling = 0;
    if (parent->ID() > PLUTO_COMPOSITE) {
	for (int j=0;j<distribution_position;j++) { //loop over all valid distributions
	    if (dist[j]->GetVersionFlag() & VERSION_MASS_SAMPLING) {
		if (num_mass_sampling)
		    dist[j]->SetEnable(0);
		num_mass_sampling++;
	    }
	}
    }
    // end check tentative distributions
    init_done=kTRUE;
    print_tentative=kFALSE;

    for (int j=0;j<distribution_position;j++) { //loop over all valid distributions
	if ((dist[j]->GetStatus() == 0) && dist[j]->GetEnable()){ //all particles set
	    dist[j]->GetDepth();  //Init sub-models
	    
	}
    }



    return kTRUE;
}

// Bool_t PChannel::InitEnv() {
//   //Init "envelope" afterburner for distributions spanning over a chain
//   cout << "**********init env" << endl; 
//   for (Int_t k=1;k<=n;k++) {
//     //looping over daughters
//     int kk=0;
//     ptcls[k]->Print();
//     while (ptcls[k]->GetDaughter(kk) != NULL) {  
//       cout << "found:" << endl;
//       ptcls[k]->GetDaughter(kk)->Print();
//       for (int j=0;j<distribution_position;j++) {
// 	dist[j]->Print();
// 	cout << dist[j]->GetStatus() << ":" << dist[j]->GetEnable() << endl;
// 	if ((dist[j]->GetStatus() == -1) && dist[j]->GetEnable()){ //not all particles set
// 	  dist[j]->SetParticle(ptcls[k]->GetDaughter(kk), 
// 			       ptcls[k]->GetDaughter(kk)->ID(), 
// 			       PARTICLE_LIST_GRANDDAUGHTER);
// 	}
//       }
//       kk++;
//     }
//   }
// }


int PChannel::Decay() {
    // The N-body phase-space decay function is based on the CERNLIB
    // routine GENBOD. 
    //
    // Return code:
    //
    // status = 0  all ok
    //        = 1  not used
    //        = 2  weighting negative
    //        = 3  not enough energy
    //        = 4  not enough energy to decay to final state
    //        = 5  genbod error (energy violation)
    //        = 6  genbod error (kinematics failed)
    //        = 7  genbod error (energy violation)
    //        = 8  EOF on input in PFileInput::readParticle()
    //        = 9  not used
    //        = 10 Abort distribution rejection after 10 times, do chain again
    //        = 11 Fermi sampling failed for deuteron wave function
    // > 80 is for PReaction
    //        = 80 Empty event
    //   STATUS_NOT_DECAYED 99

    TVector3 vCreation, vDecay;
    Double_t tCreation, tDecay;

    //    cout << "enter " << decay_key << endl;

    if (!init_done) Init();

    if (parent->IsActive() == kFALSE) {
	for (Int_t k=1;k<=n;k++) ptcls[k]->SetInActive();  // inactivate children
	return STATUS_OK;
    }
    else {
	for (Int_t k=1;k<=n;k++) ptcls[k]->SetActive();   
	// re-activate children as they are de-activated in PReaction
      
	//Vertex calculations have to be done in each decay
	vCreation = parent->GetVertex();  // creation vertex of parent
	tCreation = parent->T();          // creation time of parent
	tDecay = tCreation + parent->Gamma() * parent->GetProperTime();
	vDecay = vCreation + parent->Gamma() * parent->GetProperTime()
	    * parent->BoostVector();  // decay vertex of parent
	for (Int_t k=1;k<=n;k++) {
	    ptcls[k]->SetProperTime();           // proper decay time (in mm/c)
	    ptcls[k]->SetVertex(vDecay,tDecay);  // set creation vertex of product
	    // ptcls[k]->clearDebugString();  BUGBUG
	}
    }


  
    w=parent->W();
    int i, nProd=n;

    status=0;
    num_not_finalized = 0;

    if (thSrc) {                   // parent is a midrapidity thermal source
	ThermalSampling();
	//cout << "exit ThermalSampling" << endl;
	return status;
    }    

    if (tcSrc) {                   // parent is transport code (file input)
	status = ReadFileInput();
	return status;
    }    

    if (dlSrc) {                   // parent is a Dilepton
	MakeDilepton();
	return status;
    }

    ecm=parent->M();

    if (ecm<emin) {
	// cout << "ecm: " << ecm << " emin: " << emin <<  endl;
 	// parent->Print();
 	// Print();
	return status=4; // not enough energy to do the channel
    }

    // if (nProd==1 && ipid[0] > 1000) { // if first step is to a single resonance
    // 	if ( ecm > makeStaticData()->GetDecayEmin(ipid[1]) ) {   // ecm larger than M-5*fwhm?
    // 	    ptcls[1]->SetPxPyPzE(0.,0.,0.,ecm);    // mass of resonance = ecm
    // 	    ptcls[1]->SetW(w);
    // 	}
    // 	else  {
    // 	    return status=3;
    // 	}
    // } //removed, 22.10.2010, because I want to do this with models (IF)

    {
      
	for (int j=0;j<distribution_position;j++) { //loop over all valid distributions
	    if ((dist[j]->GetStatus() == 0) && dist[j]->GetEnable()){ //all particles set
		dist[j]->Prepare();
	    }
	}
	//	cout << "make genbod" << endl;
	Int_t do_flag=1, distr_status = 0;
	while (do_flag) {  // do genbod loop until all distributions are convinced

	    status=Genbod(nProd);  // induce first isotropic decay
	    //    	    cout << "genbod done" << endl;
	    distr_status = 0;

	    for (int j=0;j<distribution_position;j++) { //loop over all valid distributions
		if ((dist[j]->GetStatus() == 0) && dist[j]->GetEnable()){ //all particles set
		    if (dist[j]->IsValid() == kFALSE) { 
			distr_status++;
			if (dist[j]->CheckAbort()) {
			    do_flag = 1000000;//set to big number, next step will abort
			}
		    }
		  
		} //end if ((dist[j]->GetStatus() == 0) && dist[j]->GetEnable())
	    } // end for (int j=0;j<distribution_position;j++)

	    if (distr_status) {// did not work out...
		do_flag++;
		if (do_flag > 10) // abort after 10 times, do chain again
		{
		    status = 10;
		    do_flag=0;
		}
	    } else {
		do_flag=0;
	    }
	  
	} //end while (do_flag)

	for (int j=0;j<distribution_position;j++) { //loop over all valid distributions
	    if ((dist[j]->GetStatus() == 0) && dist[j]->GetEnable()){ //all particles set
	      if (!dist[j]->Finalize()) {
		//distribution is waiting for more....
		if (num_not_finalized < MAX_NUM_NOT_FINALIZED) {
		  distribution_not_finalized[num_not_finalized] = dist[j];
		  num_not_finalized++;
		} else {
		  Error("Decay","MAX_NUM_NOT_FINALIZED reached");
		}
	      }
	    }
	}
      
    }

    if (status) return status;                  // genbod failed

    //In a first step collect all weights from all distributions
    //This we have to do because in the increment we have to 
    //use the stat. scaling with 1/w
 
    Double_t new_weight=ptcls[0]->W();
    Double_t hidden_generator_weight=ptcls[0]->GenW();
    Double_t inv_generator_weight=ptcls[0]->InvGenW();
    
    Double_t dynamic_range=1.;

    if (*weight_version) {
	
#if 1
	
	Double_t generator_weight=hidden_generator_weight; 
	for (int j=0;j<distribution_position;j++) { //loop over all valid distributions
	    if ((dist[j]->GetStatus() == 0) && dist[j]->GetEnable()){ //all particles set
		dist_weight[j] = dist[j]->GetWeight();
		if ((dist[j]->GetVersionFlag()) & VERSION_INVERT_WEIGHTING) { //Here I assume a "generator"
		    generator_weight*=dist_weight[j];
		    dynamic_range*=dist[j]->GetDynamicRange();
		}
	    }
	}
	
	for (int j=0;j<distribution_position;j++) { //loop over all valid distributions
	    if ((dist[j]->GetStatus() == 0) && dist[j]->GetEnable()){ //all particles set
		if (((dist[j]->GetVersionFlag()) & VERSION_INVERT_WEIGHTING) ||
		    ((dist[j]->GetVersionFlag()) & VERSION_WEIGHTING)) {
		    
		    Double_t local_weight = dist_weight[j];
		    //weighting is enabled
		   //  dist[j]->Print();
// 		    cout << "...has " << local_weight << endl;
		
		    if ((dist[j]->GetVersionFlag()) & VERSION_GENERATOR_MC) { 
			//Pure function models
			//--> These are function for which we have
			//Int dx dsigma/dx = sigma
			//In this case, the boundaries have to be taken into account
			//This is the Delta_x in the MC integral
			//cout << "We got: " << local_weight << " with: " << dynamic_range << endl;
			local_weight*=dynamic_range;
			dynamic_range=1.;
		    }	

		    if (dist[j]->GetExpectedWeightMean() > 0) { //re-scaling is allowed
			if ((dist[j]->GetVersionFlag()) & VERSION_INVERT_WEIGHTING) {
			    local_weight=(1./local_weight);
			    dist[j]->IncrementWeightSum(local_weight);
			} else {
			    //Weighting model
			    dist[j]->IncrementWeightSum(local_weight,1./generator_weight);
			}
			//cout << dist[j]->GetExpectedWeightMean() << ":" << dist[j]->CalcWeightMean() << endl;
			local_weight*=(dist[j]->GetExpectedWeightMean()/dist[j]->CalcWeightMean());
		    } 

		    //		    cout << "1:" <<  local_weight << endl; 		   		    
		    new_weight*=local_weight;  //fold with local weights
		    //		    cout << "2:" <<  new_weight << endl; 
		    
		    dist_weight_sum[j]+=local_weight;
		    dist_counter[j]++;
		} else { //for sampling models use just the mean
		    new_weight*=dist[j]->GetExpectedWeightMean();
// 		cout << new_weight << endl;
// 		dist[j]->Print();
		    dist_weight_sum[j]+=dist[j]->GetExpectedWeightMean();
		    dist_counter[j]++;
		    
		    //If we have broad resonances, like in the case
		    //Delta->dilepton + gamma we have to keep the information
		    //about the generator
		    
		    //this is only possible if we "push" the history of used generators
		    //to the next PChannel(s)
		    
		    //check if weight is what we expect
		    Double_t local_weight = dist_weight[j];
//		hidden_generator_weight *=local_weight;
		    
		    Double_t inv_local_weight=0.;
		    if (local_weight>0.) inv_local_weight=(1./local_weight);
		    dist[j]->IncrementWeightSum(inv_local_weight);
//		inv_generator_weight*=(inv_local_weight/dist[j]->CalcWeightMean());
		    
		}
	    }
	}

	if (new_weight<=0.) return status=2; // something wrong

	weight_sum+=new_weight;
	//ptcls[0]->Print();
	//cout << "new_weight:" << new_weight << endl;
	
#endif

    }

//     for (i=1;i<=nProd;++i) {
// 	ptcls[i]->Print();
//     }

    TVector3 beta=parent->BoostVector();        // inv. mass velocity vector
    if (beta.Mag()>0.) for (i=1;i<=nProd;++i) ptcls[i]->Boost(beta);

    for (i=1;i<=nProd;++i) {
	ptcls[i]->SetParent(parent); // set pointer to parent
	//ptcls[i]->Print();
	ptcls[i]->SetW(new_weight);
	ptcls[i]->SetGenW(hidden_generator_weight);
	ptcls[i]->SetInvGenW(inv_generator_weight);

    }


    
#if 0
    for (int j=0;j<distribution_position;j++) { //loop over all valid distributions
	if ((dist[j]->GetStatus() == 0) && dist[j]->GetEnable()){ //all particles set
	    //dist[j]->WriteDebugInfo(parent);
	}
    }
#endif

    return status;
}


int PChannel:: Genbod(int nProd) {
    // N-body phase-space decay of a particle in its rest frame,
    // wrapper function to all the distribution methods

    int i; //, nm2=nProd-2;
    double conserve_e=0., em[nProd];
    Bool_t PDistribution_sampleMass = kFALSE;
    int PDistribution_sampleMomentum=0;

    //if (nm2<0) return status=5; //removed (see above)

    //MASS SAMPLING
    for (int j=0;j<distribution_position;j++) { //loop over all valid distributions
	if ((dist[j]->GetStatus() == 0) && dist[j]->GetEnable()){ //all particles set
	    //	    cout << dist[j]->GetVersionFlag() << ":" <<dist[j]->GetIdentifier() <<endl;
	    if ((dist[j]->GetVersionFlag()) & VERSION_SAMPLING)
		//		cout << dist[j]->GetIdentifier() << endl;
		if (dist[j]->SampleMass()) {
		    PDistribution_sampleMass = kTRUE;
		    for (i=0;i<nProd;++i) {
			em[i]=ptcls[i+1]->M(); //copy particle mass to array
		    }
		}
	}
    }

    conserve_e=ptcls[0]->M(); //take into account possible reconstruct

    if (PDistribution_sampleMass == kFALSE) { //No user-defined distribution    
	Warning("Genbod","No mass sampling model(s) found in %s",
		makeStaticData()->GetDecayNameByKey(decay_key));

	return status=6;
    } else { //Recalculate total energy (for user-defined mass sampling)
	for (i=0;i<nProd;++i) conserve_e-=em[i];
    }
  
    if ((conserve_e <  -1.e-8) && (nProd>1) ) {
      // cout << "***********" << endl;
      // for (i=0;i<=nProd;++i) ptcls[i]->Print();
      // cout << conserve_e << endl;
	return status=5;
    }
    conserve_e=ptcls[0]->M();

//     cout << ecm << endl;
//     ptcls[0]->Print();

    Int_t last_sampling_model = 0;

    //MOMENTUM SAMPLING
    for (int j=0;j<distribution_position;j++) { //loop over all valid distributions
	if ((dist[j]->GetStatus() == 0) && dist[j]->GetEnable()){ //all particles set
	    if ((dist[j]->GetVersionFlag()) & VERSION_SAMPLING)
		if (dist[j]->SampleMomentum()) {
		    if (!PDistribution_sampleMomentum) last_sampling_model = j;
		    PDistribution_sampleMomentum++;
		} 
	    if (PDistribution_sampleMomentum > 1) {
		Warning("Genbod","More than one momentum sampling model found in %s",
			makeStaticData()->GetDecayNameByKey(decay_key));
		Warning("Genbod","Model 1: [%s]",dist[last_sampling_model]->GetIdentifier());
		Warning("Genbod","Model 2: [%s]",dist[j]->GetIdentifier());
		
		return status=6;
	    }
	}
    }

    if (PDistribution_sampleMomentum == 0) {
	Warning("Genbod","No momentum sampling model found in %s",
		makeStaticData()->GetDecayNameByKey(decay_key));
	return status=6;
    }

    //ANGLE SAMPLING
    for (int j=0;j<distribution_position;j++) { //loop over all valid distributions
	if ((dist[j]->GetStatus() == 0) && dist[j]->GetEnable()){ //all particles set
	    if ((dist[j]->GetVersionFlag()) & VERSION_SAMPLING)
		dist[j]->SampleAngle();

	}
    }

    //check conservation of energy again
    for (i=0;i<nProd;++i) conserve_e-=ptcls[i+1]->E();

    if (fabs(conserve_e)>1.e-8) {
 	 // cout << "*************conserve_e" << conserve_e  << endl;
 	 // ptcls[0]->Print();
 	 // for (i=0;i<nProd;++i) ptcls[i+1]->Print();
 	 // cout << "*************" << endl;
	return status=7; // conservation of energy violated
    }


    else return status=STATUS_OK;
}




int PChannel::SetDistribution(PDistribution * distribution) {
    //return -1 if distribution does not match the channel

    if (distribution->GetEnable()==0) return -1;

    if (distribution_position == MAX_DISTRIBUTIONS) {
	Warning("SetDistribution","MAX_DISTRIBUTIONS reached");
	return -1;
    }

    Int_t check=0;

    PParticle *p0=ptcls[0];
    
    //The grand(grand)parent and sisters not known yet at this stage, 
    //will be done "online"
    //The only thing what I can do here is to have a 1st rough look

    distribution->Reset();

//    if (p0->ID() < 1000) //skip quasi-particles
//	check = distribution->SetParticle(p0, p0->ID(), PARTICLE_LIST_PARENT); 
//    else
//	check = distribution->SetParticle(p0, DISTRIBUTION_QUASI_PID, PARTICLE_LIST_PARENT); 

//either one of these methods has to work:
    check = distribution->SetParticle(p0, p0->ID(), PARTICLE_LIST_PARENT); 
    if (check<0 && (p0->ID() > 1000))
	check = distribution->SetParticle(p0, DISTRIBUTION_QUASI_PID, PARTICLE_LIST_PARENT); 
	

    for (int i=0;i<n;i++) {
	check += distribution->SetParticle(ptcls[i+1], ptcls[i+1]->ID(), PARTICLE_LIST_DAUGHTER); 
    }
   
    check += distribution->CheckDaughters(); //Are we still a daughter missing?
 
    //now check if this worked out....
    if (check < 0) {
	//forget it!
	if (distribution->debug_flag)
	     distribution->Print();

	distribution->Reset();
	return -1;
    }
    //success

    //Clone the distribution (done before in PDistributionManagerUtil)
    distribution = (PDistribution*) distribution->Clone();

    dist[distribution_position]=distribution;
    distribution_position++;
    if (distribution->GetStatus() == 0) { //we are complete
	if (distribution->Init() == kFALSE) {
	    printf("                             [error] %s:%s\n",distribution->GetDescription(),distribution->GetIdentifier());
	    distribution->SetEnable(0);
	} else distribution->GetDepth();  //Init sub-models
    } else { //if we are not yet complete, relatives are used
	if (init_done) { //relatives are already available
	    Reset();
	    Init();
	    if (distribution->GetStatus() == 0) { //we are complete
		if (distribution->Init() == kFALSE) {
		    printf("                             [error] %s:%s\n",distribution->GetDescription(),distribution->GetIdentifier());
		    distribution->SetEnable(0);
		} else distribution->GetDepth();  //Init sub-models
		
		//printf("                             [incomplete] %s\n",distribution->GetDescription());
	    }
	} //END init_done
    } //...end not yet complete

       
    return 0;
}



void PChannel:: GetMessage() { printf(" %s\n",Message[status]); }

char const *PChannel::GetName(void) const {
    //returns the name from the data base
    if (!thSrc && !tcSrc && !dlSrc) {
	//Here, the string from PData should be used!!!
	if (decay_key < 0 ) {
	    Warning("GetName","No decay key found");
	    exit(1);
	}
	return makeStaticData()->GetDecayNameByKey(decay_key);
    }
    if (thSrc) return " Fireball";
    else if (tcSrc) return " File Input";
    else if (dlSrc) return " Dilepton source";
    return "<unknown>";
}

void PChannel::PrintReport() const {
    PrintReaction();
    cout << "Weighting report" << endl;
    for (int j=0;j<distribution_position;j++) { //loop over all valid distributions
	printf("        [%s] %s\n",dist[j]->GetIdentifier(),dist[j]->GetDescription());

// 	if (((dist[j]->GetVersionFlag()) & VERSION_INVERT_WEIGHTING) ||
// 	    ((dist[j]->GetVersionFlag()) & VERSION_WEIGHTING)) {

	    if (dist[j]->GetExpectedWeightMean() > 0) { //re-scaling is allowed
		cout << "Expected mean: " << dist[j]->GetExpectedWeightMean() << " Reached: " 
		     << dist_weight_sum[j]/dist_counter[j] << endl;
		cout << "Calc mean: " << dist[j]->CalcWeightMean() << endl;
	    } else {
		cout << "Mean: " << dist_weight_sum[j]/dist_counter[j] << endl;
	    }
	    
//	}
    }
}

void PChannel::PrintReaction(Int_t check_key) const {

    if (check_key) {
	if (!thSrc && !tcSrc && !dlSrc) {
	    if (decay_key < 0 ) {
		Warning("PrintReaction","No decay key found");
		exit(1);
	    }
	    cout << makeStaticData()->GetDecayNameByKey(decay_key) << endl;
	    return;
	}
    }

    int j;

    if (thSrc) printf("Fireball");
    else if (tcSrc) printf("File Input");
    else if (dlSrc) printf("Dilepton source");
    else printf("%s", makeStaticData()->GetParticleName(ipid[0]));
    printf(" --> %s",makeStaticData()->GetParticleName(ipid[1]));
    for (j=2;j<=n;++j) printf(" + %s",makeStaticData()->GetParticleName(ipid[j]));
    cout << endl;

}

void PChannel::PrintNew() {
    
    Long_t nPattern=0;
    for (int j=0;j<distribution_position;j++)
	if (dist[j]->GetEnable() != 0) {
	    nPattern += ((j<64) ?  1<<(j+1) : 0);
	}
    if (nPattern & (~fEnablePattern)) { 
	fEnablePattern |= nPattern;
	Print();
    }
}

void PChannel::Print(const Option_t* delme) const {
    int j;

    if (quasi_pchannel)
 	quasi_pchannel->Print();

    PrintReaction();
    printf("        Interaction model(s):\n");
    //Distributions
    for (j=0;j<distribution_position;j++) {
	
	if (dist[j]->GetEnable() != 0) {
	    if (dist[j]->GetStatus() == 0) { //all particles set 
		dist[j]->GetDepth();  //Init sub-models
		if (dist[j]->Path())
		    printf("        [%s] %s {/%s}",dist[j]->GetIdentifier(),dist[j]->GetDescription(),dist[j]->Path());
		else
		    printf("        [%s] %s",dist[j]->GetIdentifier(),dist[j]->GetDescription());
		dist[j]->SubPrint(0);
		cout << endl;
	    }
	    else
		printf("        [%s,tentative] %s\n",dist[j]->GetIdentifier(),dist[j]->GetDescription());
	    if (delme != NULL && strcmp(delme,"debug") == 0) {
		dist[j]->Print();
	    }
	    
	}
    }
    if (thSrc) ((PFireball * )parent)->Print("      ");
}

ClassImp(PChannel)
