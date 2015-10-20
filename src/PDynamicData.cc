/////////////////////////////////////////////////////////////////////
//  Pluto Dynamic Data Wrapper
//
//  Provides wrapper functions to deal with coupled channel
//  calculations
//
//
//                             Author:  IF
//                             Written: 23.7.2007
//                             Revised: 
//
//////////////////////////////////////////////////////////////////////

#include "PDynamicData.h"
#include "PStdData.h"
#include "PUtils.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "Pdefines.h"


PDynamicData& fDynamicData()
 {
   static PDynamicData* ans = new PDynamicData();
   return *ans;
 }

PDynamicData * makeDynamicData()
 {
   return &fDynamicData();
 }

PDynamicData::PDynamicData() {

    PDataBase *base=makeDataBase();
    i_result=NULL;
    c_result=NULL;
    d_result=NULL;
    t_result=NULL;
    num_pchannels=0;

    pid_param=base->GetParamInt("pid");
    name_param=base->GetParamString("name");
    model_param=base->GetParamTObj("model");
    didx_param=base->GetParamInt("didx");
    scfactor_param=base->GetParamDouble("scfactor");
    sccount_param=base->GetParamInt("sccount");

    pnmodes_param = base->GetParamInt ("pnmodes");
    link_param    = base->GetParamInt ("link");



    Info("PDynamicData()","(%s)", PRINT_CTOR);
    
  
}

Double_t PDynamicData::GetDecayPartialWidth(Double_t mass,Int_t didx) {
    //returns the partial decay width of the decay channel "didx"
    //for a given "mass"
    //If no primary channel model has been found
    //return the static width
    Double_t width;

    //Check for threshold
    if (makeStaticData()->GetDecayEmin(didx)>mass) return 0;

    PChannelModel *model = GetDecayModel(didx);
    if (!model) {
	cout << "model not found for " << didx << endl;
	return makeStaticData()->GetDecayPartialWidth(didx);
    }
    Double_t scfactor=GetDecaySCFactor(didx);
    if (model->GetWidth(mass, &width)) return (width*scfactor);
    return makeStaticData()->GetDecayPartialWidth(didx)*scfactor;
}



double PDynamicData:: GetParticleTotalWidthSum(Double_t mass,Int_t id, Int_t flag) {
    // Total decay width by pid & invariant mass (GeV/c^2).
    //
    // For most purposes GetParticleTotalWidth may be used for the total width.
    // GetParticleTotalWidthSum should be used instead only for sampling BR's.
    // Due to roundoff error in the interpolation, Width() returns a value
    // that differs slightly from the sum of the partial widths Width1().
    // As a result, small branching ratios (e.g. electromagnetic processes)
    // may never be sampled if Width() is used. For this reason, if a BR is 
    // requested, it is best to calculate the total width explicitly by 
    // summing the partial widths for any given mass.
    // 
    // flag=1: take into account only hadronic decays

    static double mo=0.;
    static double wt;
    static int io=0;
    static int old_flag=0;

    if (!makeStaticData()->IsParticleValid(id)) { // pid out of range
	cout << PRINT_WARNING << "PData:: GetDecayPartialTWidth: id " 
	     << id << " out of range" << endl;
	return 0.;
    }

    int twidx=makeStaticData()->GetTWidx(id);

    if (!twidx) GetParticleDepth(id);  // initialize indices on first call

    if (twidx==-1) {
	mo=mass;
	io=id;
	return wt=makeStaticData()->GetParticleTotalWidth(id);
	// static total width
    }
    
    if (io!=id) {                        // parent id changed since last call
	io=id;
	mo=0.;
    }

    if (mo!=mass || flag!=old_flag) {     // mass changed since last call
	mo=mass;
	old_flag=flag;
	if (mass<PData::LMass(id) || mass>PData::UMass(id)) 
	    return 0.;    // out of range

	double g_tmp=0.;         // reset sums

	//now loop over decay modes
	Int_t key=makeDataBase()->GetEntryInt("pid", id);
	Int_t *didx;
	Int_t listkey=-1;
	while (makeDataBase()->MakeListIterator
	       (key, "pnmodes", "link", &listkey)) {
	    makeDataBase()->GetParamInt 
		(listkey, "didx" , &didx); //small workaround
	    if (makeStaticData()->GetPWidx(*didx)==-1
		&&mass>=makeStaticData()->GetDecayEmin(*didx)) 
		// if mode is unknown but kinematically accessible...
		g_tmp+=
		    makeStaticData()->GetDecayBR(*didx)*
		    makeStaticData()->GetParticleTotalWidth(id);
	    // ...update with the current static BR
	    else if (flag==0) 
		g_tmp+=GetDecayPartialWidth(mass,*didx);         
	    else if (makeStaticData()->IsDecayHadronic(*didx))
		g_tmp+=GetDecayPartialWidth(mass,*didx);  
	    // otherwise, decay width of known mode
	}
	wt=g_tmp;          
    }

    return wt;
}


double PDynamicData:: GetDecayBR(int idx, double m) {
    // returns branching ratio by mode index and mass (GeV/c^2)

    PChannelModel *model = GetDecayModel(idx);
    if (!model) {	
	return makeStaticData()->GetDecayBR(idx);
    }
    Double_t br;
    Double_t tw=GetParticleTotalWidth(m,makeStaticData()->GetDecayParent(idx));
    if (model->GetBR(m, &br, tw)) return br;
    return makeStaticData()->GetDecayBR(idx);

// 	return (wt>0.) ? wp/wt : 0.;
//->return (wt>0.) ? PBR[idx]*PWidth[id]/wt : 0.; <- alternative
}


void PDynamicData::ListDecayBR(int id, double m) {
    // list branching ratios of particle id for mass m
    
    if (!makeStaticData()->IsParticleValid(id)) {  // pid out of range
	return;
    }
    printf("\n%s branchings for mass=%f GeV   [low|pole|up|w|w0]=[%f|%f|%f|%e|%e] GeV\n\n",
	   makeStaticData()->GetParticleName(id),PData::LMass(id),
	   makeStaticData()->GetParticleMass(id),
	   makeStaticData()->GetParticleMass(id),PData::UMass(id),
	   makeStaticData()->GetParticleTotalWidth(id),GetParticleTotalWidthSum(m,id));
    if (m<PData::LMass(id) || m>PData::UMass(id)) { // mass out of range
	Warning("ListDecayBR","mass out of range");
	return;
    }

    //now loop over decay modes
    Int_t key=makeDataBase()->GetEntryInt("pid", id);
    Int_t listkey=-1;
    int i=1;

    double sum=0;

    while (makeDataBase()->MakeListIterator
	   (key, "pnmodes", "link", &listkey)) {
	int pos=makeStaticData()->GetDecayIdxByKey(listkey);
	double br=GetDecayBR(pos,m);  // compare with BR of current mode
	sum+=br;
	double pw = GetDecayPartialWidth(m,pos);
	printf("%d: idx=%d sw=%2d  width=%e GeV    br=%f %%  (Width1()=%e)\n",
	       i,pos,makeStaticData()->GetPWidx(pos),
	       makeStaticData()->GetDecayPartialWidth(pos),
	       100.*br, pw);
	i++;
    }
    printf("\nTotal width=%e GeV   Sum(br)=%f %%\n",GetParticleTotalWidthSum(m,id),100.*sum);

}

int PDynamicData:: PickDecayChannel(const int & id, const double & m) {
    // returns the index of a decay mode for particle pid=id of mass m (GeV/c**2)
    // selected randomly, consistent with the branching ratios

    
    if (!makeStaticData()->IsParticleValid(id)) {
	// pid out of range
	cout << PRINT_WARNING << "PDynamicData:: PickDecayChannel: id " << id << " out of range" << endl;
	return -1;
    } else if (m<PData::LMass(id) || m>PData::UMass(id)) return -2; // mass out of range
    
    int n=makeStaticData()->
	GetParticleNChannels(id); // number of channels available
    
    if (!n) return -3;                     // stable particle
    double r=PUtils::sampleFlat(), br=0.;  // pick a random number
    double sum=0.;
    Int_t key=makeDataBase()->GetEntryInt(pid_param, id);
    Int_t *didx;
    Int_t listkey=-1;

    while (makeDataBase()->MakeListIterator
	   (key, pnmodes_param, link_param, &listkey)) {	
	makeDataBase()->GetParamInt 
	    (listkey, didx_param , &didx); //small workaround -> should work on keys
	sum+=GetDecayBR(*didx,m); // normalize first branching ratios
	//cout << "BR is " << GetDecayBR(*didx,m) << " for " << didx[0] << endl;
    }

    if (sum<=0.) {
	Warning("PickDecayChannel","sum=%f\n",sum);
	return -4;
    }
    listkey=-1;
    while (makeDataBase()->MakeListIterator
	   (key, pnmodes_param, link_param,  &listkey)) {
	if (!makeDataBase()->GetParamInt 
	    (listkey, didx_param , &didx)) {
	    Warning("PickDecayChannel","didx not found for listkey %i",listkey);
	}
	br+=GetDecayBR(*didx,m);    // normalize first branching ratios
	if (r<br/sum) {
	    return *didx; // return selected index
	}
    }
    
    Warning("PickDecayChannel","id=%d sum=%f\n",id,sum);
    return -5;
}



int PDynamicData:: PickDecayChannel(const int & id, const double & m, int * array) {
    // returns the number and pids of the decay products for the decay 
    // of particle pid=id of mass m (GeV/c**2), via a mode selected
    // randomly, consistent with the branching ratios
    // In addition, it does not hurt to return the index (IF)

    if (!array) {  // invalid address of array
	Warning("PickDecayChannel","Invalid address of array");
	return -1;
    } else if (!makeStaticData()->IsParticleValid(id)) {            
	// pid out of range
	Warning("PickDecayChannel","id %i out of range",id);
	return -1;
    }
    
    array[0]=0;  // failed to select an index
    Int_t p = PickDecayChannel(id,m);
    array[0]=10; //max 10 particles -> for the Mode() BUGBUG: Size of array not checked
    if (p>=0) makeStaticData()->GetDecayMode(p,array); 
    // return number of products and pids
    return p;
}

Double_t PDynamicData::GetParticleScalingFactor(Int_t didx) {
    if (! makeDataBase()->GetParamDouble (pid_param, didx , scfactor_param, &d_result)) {
	return 1.;
    }
    return *d_result;
}

void PDynamicData::SetParticleScalingFactor(Int_t didx, Double_t factor) {
    if (! makeDataBase()->GetParamDouble (pid_param, didx , scfactor_param, &d_result)) {
	makeDataBase()->SetParamDouble (makeDataBase()->GetEntryInt(pid_param, didx) , "scfactor", new Double_t(factor));
	return;
    }
    *d_result *=factor;
}


Double_t PDynamicData::GetDecaySCFactor(Int_t didx) {
    if (! makeDataBase()->GetParamDouble (didx_param, didx , scfactor_param, &d_result)) {
	return 1.;
    }
    return *d_result;
}

void PDynamicData::SetDecaySCFactor(Int_t didx, Double_t factor) {
    if (! makeDataBase()->GetParamDouble (didx_param, didx , scfactor_param, &d_result)) {
	cout << PRINT_WARNING << "PDynamicData::SetDecaySCFactor: factor not present" << endl;
	return;
    }
    *d_result *=factor;
}

bool PDynamicData::CheckSCAbort(Int_t didx) {
    //to stop endless loops
    if (! makeDataBase()->GetParamInt (didx_param, didx , sccount_param, &i_result)) {
	makeDataBase()->SetParamInt (makeDataBase()->GetEntryInt(didx_param, didx) , "sccount", new Int_t(10)); //TODO:MAX_TRIES
	return kTRUE;
    }
    (*i_result)--;
    if (*i_result <= 0) return kFALSE;
    return kTRUE;
}

PChannelModel *PDynamicData::GetDecayModel(Int_t didx) {
    //returns the primary decay model

    if (! makeDataBase()->GetParamTObj (didx_param, didx , model_param, &t_result)) {
	return NULL;
    }
    return (PChannelModel *) t_result;
}

PChannelModel *PDynamicData::GetDecayModelByKey(Int_t key) {
    //returns the primary decay model

    if (! makeDataBase()->GetParamTObj (key , model_param, &t_result)) {
	return NULL;
    }
    return (PChannelModel *) t_result;
}

PChannelModel *PDynamicData::GetDecayModelByKey(Int_t key, Int_t defkey) {
    //returns the secondary decay model
    

    Int_t listkey = makeStaticData()->GetSecondaryKey(key,defkey);
    if (listkey<0) return NULL;

    if (makeDataBase()->GetParamTObj (listkey , model_param, &t_result)) {
//	    cout << (((PChannelModel *) t_result) -> GetDef()) << endl;
	if ((((PChannelModel *) t_result) -> GetDef())  == defkey)
	    return (PChannelModel *) t_result;
	else {
	    Warning("GetDecayModelByKey","Consistency check failed");
	    return NULL;
	}
    }
    return NULL;

}



PChannelModel *PDynamicData::GetParticleModel(Int_t pid) {
    //returns the primary particle model

    if (! makeDataBase()->GetParamTObj (pid_param, pid , model_param, &t_result)) {
	return NULL;
    }
    return (PChannelModel *) t_result;
}


PChannelModel * PDynamicData::GetParticleSecondaryModel(const char * name, 
							const char * modelname) {
    Int_t sec_key = makeStaticData()->MakeDirectoryEntry("modeldef",NMODEL_NAME,LMODEL_NAME,modelname);
    if (sec_key<0) return NULL;

    Int_t primary_key = makeStaticData()->GetParticleKey(name);
    Int_t listkey = makeStaticData()->GetSecondaryKey(primary_key,sec_key);
    if (listkey<0) return NULL;

    Int_t model_param = makeDataBase()->GetParamTObj("model");
    TObject *t_result = NULL;
    if (makeDataBase()->GetParamTObj (listkey , model_param, &t_result)) {
	if ((((PChannelModel *) t_result) -> GetDef())  == sec_key)
	    return (PChannelModel *) t_result;
	else {
	    Warning("GetSecondaryModel","Consistency check failed");
	    return NULL;
	}
    }
    return NULL;
}

Double_t PDynamicData::GetParticleTotalWidth(Double_t mass,Int_t pid) {
    //returns the particle width if a models exists, if not the static
    //width
    PChannelModel *model=GetParticleModel(pid);
    Double_t width;
    if (model) {
	if (model->GetWidth(mass, &width))
	    return width;
    }
    return makeStaticData()->GetParticleTotalWidth(pid);
}

Double_t PDynamicData::GetParticleTotalWeight(Double_t mass,Int_t pid, Int_t didx) {
    //return the mass-dependent width if a model exists
    //no mode: use the fixed Breit-Wigner with the total width
    
    PChannelModel *model=GetParticleModel(pid);

    if (model) {
	Double_t w = model->GetWeight(&mass,&didx);	
	if (w>0) return w;
	return 0;
    }
    
    return 0; //BUGBUG make fixed BW here

}


bool PDynamicData::SetDecayModel(Int_t didx, PChannelModel * model) {
    Int_t key=makeDataBase()->GetEntryInt("didx", didx);
    if (! makeDataBase()->SetParamTObj (key,"model", (TObject *)model)) {
	return kFALSE;
    }   
    return kTRUE;
}

bool PDynamicData::SetDecayModelByKey(Int_t didx, PChannelModel * model) {

    if (! makeDataBase()->SetParamTObj (didx,"model", (TObject *)model)) {
	return kFALSE;
    }   
    return kTRUE;
}


double PDynamicData:: GetParticleLife(const int & id, double m, int idx) {
    // mean life
    // Arguments: 1. id=particle id 
    //            2. m=mass (GeV/c**2)
    //            3. idx=decay-mode index (PPosition)
    
    //CALLED FROM: PParticle, PReaction
    const long double hbar=6.582122e-25;// (s GeV/c**2)

    double w0 = 
	makeStaticData()->GetParticleTotalWidth(
	    makeStaticData()->IsParticleValid(id)); // static total width
    if (id==50 || 
	makeStaticData()->IsParticle(id,"dilepton")) 
    {  // dilepton = virtual particle!
	if (w0 > 0.0) return hbar/w0;
	else if (m > 0.0) return hbar/m;
	else return 0.0;
    }

    double tau0=(w0>0.) ? hbar/w0 : 1.e16; 
    // static total mean life (width=0 means stable!)
  
    if (m>0.&&w0>0.) {  // unstable particle with non-zero mass
	if (idx==-1) {  // request for total mass-dependent mean life
	    double w=GetParticleTotalWidth(m,id);     // total mass-dependent width
	    if (w>0.) return hbar/w;  
            // total mass-dependent mean life, if non zero width
	} else {     // request for partial mass-dependent mean life
	    double w1=GetDecayPartialWidth(m,idx);     
	    // partial mass-dependent width for decay via channel idx
	    if (w1>0.) return hbar/w1;      
	    // partial mass-dependent mean life, if non zero partial width
	}
    }
    return tau0;   // if none of the above return the static mean life
}

int PDynamicData::GetParticleDepth(const int & id, int flag) {
    // Returns the total number of embedded recursive decays 
    // (TDepth if flag=0) or the number of embedded recursive
    // hadronic decays (HDepth if flag=1).
    // The default value 0 means the index has not been accessed yet.
    // Zero depth corresponds to index -1.
  
    if (makeStaticData()->IsParticleValid(id) == 0) {
	cout << "PDynamicData::GetParticleDepth: Particle id out of range " << id << endl;
	return 0;
    } else if (makeStaticData()->GetTDepth(id)) {
	return (flag) ? (
	    (makeStaticData()->GetHDepth(id)!=-1) ? 
	    makeStaticData()->GetHDepth(id) : 0) :
	    ((makeStaticData()->GetTDepth(id)!=-1) ? 
	     makeStaticData()->GetTDepth(id) : 0);
    }
    // Enters once on the first call;
    // It also turns off the full-width (TWidx) and partial-width
    // indices (PWidx) for particles and decay modes, respectively,
    // if the total and partial widths are not calculated. In the 
    // latter case the static width and branching ratio are used.

    int n=makeStaticData()->GetParticleNChannels(id);
    
    if (!n) {                             // no decay modes
	
	makeStaticData()->SetTWidx(id,-1);       // turn off total-width flag
	makeStaticData()->SetTDepth(id,-1);      // turn off total depth index
	makeStaticData()->SetHDepth(id,-1);      // turn off hadronic depth index
	
    } else { 
	int iter=0, iter2=0,  np;
	
	//now loop over decay modes
	Int_t key=makeDataBase()->GetEntryInt("pid", id);
	Int_t listkey=-1;
	
	Int_t tid[11];
	while (makeDataBase()->MakeListIterator
	       (key, "pnmodes", "link", &listkey)) {
	    tid[0]=10; 
	    makeStaticData()->GetDecayModeByKey(listkey,tid); // retrieve current mode info
	    int pos=makeStaticData()->GetDecayIdxByKey(listkey);
	    np=tid[0];

	    //Get primary model, --> question shifted to local
	    PChannelModel *model = GetDecayModel(pos);
	    if (!model) {
		cout << "PDynamicData::GetParticleDepth: pwidx -1" << endl;
		makeStaticData()->SetPWidx(pos,-1); //switch off
	    } else {
		iter=TMath::Max(iter,   model->GetDepth()+1);
		iter2=TMath::Max(iter2, model->GetDepth(1)+1);
	    }

	}


	if (!iter&&!iter2)
	    makeStaticData()->SetTWidx(id,-1);// do not compute total width
	makeStaticData()->SetTDepth(id,(iter)?iter:-1);          // store total depth
	makeStaticData()->SetHDepth(id,(iter2)?iter2:-1);        // store hadronic depth
    }
    return (flag) ? 
	((makeStaticData()->GetHDepth(id)!=-1) ? 
	 makeStaticData()->GetHDepth(id) : 0) :
	((makeStaticData()->GetTDepth(id)!=-1) ? 
	 makeStaticData()->GetTDepth(id) : 0);
}


void  PDynamicData::PrintParticle(int pid) {
    return PrintParticleByKey(makeDataBase()->GetEntryInt("pid",pid));
};


void  PDynamicData::PrintParticleByKey(int key) {
    makeDataBase()->ListEntries(key,0,"name,pid");

    if (!makeStaticData()->GetParticleNChannelsByKey(key)) {
	cout << "This particle is stable" << endl;
	return;
    }

    PChannelModel *m = GetParticleModel(makeStaticData()->GetParticleIDByKey(key));

    if (m) m->Print();
    else cout << "<No particle model defined>" << endl;

    cout << "This particle decays via the following modes:" << endl;

    //now loop over decay modes
    Int_t listkey=-1;

    while (makeDataBase()->MakeListIterator(key, "pnmodes", "link", &listkey)) {
	PrintDecayByKey(listkey);
	m=GetDecayModelByKey(listkey);
	if (m) m->Print();
	    //m->Dump();
    }
  
};

void  PDynamicData::PrintDecayByKey(int key) {
    makeDataBase()->ListEntries(key,0,"name,didx");
    
};


PParticle *PDynamicData::GetBatchParticle(const char * name, Int_t make_val) {
    Int_t key_a = makeStaticData()->MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,name);

    TObject *ret;
    Int_t batch_particle_param = makeDataBase()->GetParamTObj("batch_particle");
    if (batch_particle_param<0) 
	batch_particle_param = makeDataBase()->MakeParamTObj("batch_particle", "PParticle storage for batch");

    if (!makeDataBase()->GetParamTObj(key_a ,batch_particle_param,&ret)) {
	if (make_val) {
	    PParticle * delme =  new PParticle("dummy",0);
	    makeDataBase()->SetParamTObj(key_a,"batch_particle", delme);
	    return delme;
	} else return NULL;
    }
    return (PParticle *)ret;
};

ClassImp(PDynamicData)
