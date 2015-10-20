////////////////////////////////////////////////////////
//  Pluto bulk decay base 
//
//  This class let all particles with tau<tauMax 
//  decay in the pluto way
//
//                    Author:  Ingo Froehlich
//                    Written: 10/07/2007
//                    Revised:
//
////////////////////////////////////////////////////////

#include "PPlutoBulkDecay.h"
#include "PDynamicData.h"
#include "PChannel.h"
#include "PDistributionManager.h"



PPlutoBulkDecay::PPlutoBulkDecay() {
    makeStdData();
    didx=makeDataBase()->GetParamInt("didx");
    stackchannel=makeDataBase()->GetParamTObj("stackchannel");

    fPriority=DECAY_PRIORITY;
}

bool PPlutoBulkDecay::Modify(PParticle ** stack, int *decay_done, int * num, int stacksize) {
    //decay all particle
    //input: particle array with *num members
    //new particles should be instantiated and *num increased

    const int maxnp = 7; //max 7 particle decay -> Matches to static-data (IF)
    Int_t st_i1, st_i2, st_i3, k;
    Int_t size_ini=*num;  //Input size


    Int_t child_ids[maxnp+1];
    PParticle *cur_p, *work[maxnp+1];
    //cout << *num << endl;

    //I follow the LIFO principle, starting with the last particle and going to 0
    //This is needed because I first HAVE to end a decay chain, before re-using the PChannel
    //BUGBUG: This might not work with a self-decay like Delta->Delta+gamma
    st_i1 = size_ini-1;

    
    st_i2 = st_i3 = size_ini;

    while (st_i1 > -1) {
	cur_p = stack[st_i1];  //set current pointer
	// 	cout << cur_p->IsActive()  << ":"<<  makeDynamicData()->GetParticleLife(cur_p->ID() ) << endl;

	if (cur_p->IsActive() && 
	    makeDynamicData()->GetParticleLife(cur_p->ID()) 
	    < tauMax && !decay_done[st_i1]) {

// 	    cout << "Current stack particle: " << st_i1 << endl;
// 	    cur_p->Print();

	    //	    cout << "cur_p" << cur_p->ID()<< ":"<< cur_p->M()<< endl;

	    Int_t channel_idx = 
		makeDynamicData()->PickDecayChannel(cur_p->ID(), cur_p->M(), child_ids);	    

	    if (channel_idx>-1) { 
//		cout << "channel_idx" <<channel_idx<< endl;
		Int_t np = child_ids[0];   // number of decay products
		if (np > 0) { // np==0 can happen if M()<threshold
		    //First we check if we have already the PChannel on stock
		    TObject * ch;
		    PChannel* d_ch;
		    
		    if (!makeDataBase()->GetParamTObj (didx, channel_idx , stackchannel, &ch)) {
			Info("Modify","(%s) Making new channel for didx %i", PRINT_AUTO_ALLOC,channel_idx);

			if (np > maxnp) {
			    printf("PPlutoBulkDecay::Modify(): np=%d > maxnp=%d\n",np,maxnp);
			    np = maxnp;
			}

			for (k=1; k<=np; k++) {
			    work[k] = new PParticle(child_ids[k]);
			}
			work[0] = new PParticle(cur_p->ID());

			Int_t sourceIdSave = cur_p->GetSourceId();  // Need to save because
			Int_t parentIdSave = cur_p->GetParentId();  // PChannel constructor
			PParticle ** local= new PParticle*[np+1];
			for(int i=0; i<=np; i++) local[i] = work[i]; //make local copy as work array is overwritten

			d_ch = new PChannel(local,np,1,1,0); // overwrites values.

			makeDistributionManager()->Attach(d_ch);
			d_ch->SetPrintTentative(0);
			d_ch->Init();

			work[0]->SetSourceId(sourceIdSave);           // Reset them here.
			work[0]->SetParentId(parentIdSave);

			makeDataBase()->SetParamTObj(makeDataBase()->GetEntryInt("didx", channel_idx) , 
						     "stackchannel", d_ch);	
			
		    }

		    if (makeDataBase()->GetParamTObj (didx, channel_idx , stackchannel, &ch)) {
			d_ch=(PChannel*)ch;
			d_ch->Reset(); //make channel clean for next use

			(d_ch->GetParticles()[0])=cur_p;

			work[0]=cur_p;

			//Get out the daughters and modify them
			for (k=1; k<=np; k++) {
			    work[k] = (d_ch->GetParticles()[k]);  //re-use daughter particles
			    work[k]->SetActive();
			    work[k]->SetIndex(st_i3+1+k);
			    work[k]->SetParent(work[0]);
			    decay_done[st_i3+k-1]=0;
			}
			for (k=1; k<=np; k++) {
			    //set the siblings again
			    if (k!=np) {
				work[k]->SetSibling(work[k+1]);
			    }
			    else {
				work[k]->SetSibling(work[1]);
			    }
			}
			
			//--> otherwise we are not doing sampling with the total width
		    } 
		    else {
			Fatal("Modify","Unable to get PChannel" );
			return kFALSE;
		    }

#if 0
		    //before doing the decay, maybe dump the stack
		    cout << "dumping stack" << endl;
		    for (int t=0;t<st_i3+np;t++) {
			cout << t << ":" << stack[t] << ":"<< stack[t]->GetSibling() << endl;
			stack[t]->Print();
		    }
#endif
		    //end debug
		    

		    (d_ch->GetParticles()[0])->SetParent(work[0]->GetParent());
		    //now Init in background
		    d_ch->SetPrintTentative(0);
 		    d_ch->Init();

		    d_ch->PrintNew();
		    d_ch->SetPrintTentative(1);

		    Int_t error_count=0;
		    if (error_count<100) {
			if (d_ch->decay()) {
			    error_count++;
			}
		    } else {
			Warning("Modify","Decay not possible, giving up");
		    }

		    work[0]->SetDecayModeIndex(channel_idx,1); //set the index, but destroy it later
		    //-->has to happen after decay()

		    decay_done[st_i1]=1;  
		    //never let this particle decay again in this event

		    //after the decay we put all local particles to stack
		    for(int k=1; k<=np; k++) {
			*stack[st_i3+k-1] = *work[k];
		    }

		    for(int k=1; k<=np; k++) {
			stack[st_i3+k-1]->SetParent(cur_p);
			stack[st_i3+k-1]->SetActive();

			//set the siblings again
			if (k!=np) {
			    stack[st_i3+k-1]->SetSibling(stack[st_i3+k]);
			}
			else {
			    stack[st_i3+k-1]->SetSibling(stack[st_i3]);
			}
		    }
		    st_i3+=np;
		    if (recursiveMode)
			st_i1=st_i3;
		    else
			st_i1=size_ini-1;

		} //if (np > 0)
	    } //idx>-1

	    if(st_i3>(stacksize-maxnp))     // check for stack overflow, taking into account next particles
		printf("PPlutoBulkDecay::decayAll: st_i3=%d > stacksize\n",st_i3);
	
	} // if (cur_p->isActive()...

	//set the focus to the last particle, and let us drop to the first one
	st_i1--;
    }

    *num=st_i3;  //Output size
	    
    return kTRUE;
}

ClassImp(PPlutoBulkDecay) 
