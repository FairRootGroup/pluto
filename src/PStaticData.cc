/////////////////////////////////////////////////////////////////////
//  Pluto Data Wrapper
//  Provides a lot of helpfull wrapper function to deal with the data
//  base. 
//
//
//                             Author:  M.A. Kagarlis
//                             Written: 31/01/99
//                             Revised: 24/05/2004 R.H.
//                             Reimplemented: 3/5/2007 IF
//
//////////////////////////////////////////////////////////////////////

#include "PStaticData.h"
#include "PStdData.h"
#include "PUtils.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "Pdefines.h"

extern char*version_string;

PStaticData& fStaticData()
{
    static PStaticData* ans = new PStaticData();
    return *ans;
}

PStaticData * makeStaticData()
{
    return &fStaticData();
}

PStaticData::PStaticData() {
    if (!makeStdData()->fillDataBase()) {
	Error ("PStaticData","Data base could not be filled");
	exit(1);
    }
    PDataBase *base=makeDataBase();
    i_result=NULL;
    c_result=NULL;
    d_result=NULL;
    t_result=NULL;
    
    freeze=kFALSE;

    pid_param    = base->GetParamInt("pid");
    name_param   = base->GetParamString("name");
    meson_param  = base->GetParamInt("meson");
    baryon_param = base->GetParamInt("baryon");
    lepton_param = base->GetParamInt("lepton");
    charge_param = base->GetParamInt("charge");
    spin_param   = base->GetParamInt("spin");
    ispin_param  = base->GetParamInt("ispin");
    parity_param = base->GetParamInt("parity");
    mass_param   = base->GetParamDouble("mass");
    width_param  = base->GetParamDouble("width");
    pkf_param    = base->GetParamInt("pkf");
    didx_param   = base->GetParamInt("didx");
    widx_param   = base->GetParamInt("widx");
    mesh_param   = base->GetParamTObj("mesh");
    tf1_param    = base->GetParamTObj("tf1");
    ethreshold_param = base->GetParamDouble("ethreshold");
    lmass_param  = base->GetParamDouble("lmass");
    umass_param  = base->GetParamDouble("umass");
    brflag_param = base->GetParamInt("brflag");
    tdepth_param = base->GetParamInt("tdepth");
    hdepth_param = base->GetParamInt("hdepth");
    br_param     = base->GetParamDouble("br");
    brorig_param = base->GetParamDouble("brorig");
    count_param  = base->GetParamInt("pnmodes");
    d1_param     = base->GetParamInt("d1");
    d2_param     = base->GetParamInt("d2");
    d3_param     = base->GetParamInt("d3");
    d4_param     = base->GetParamInt("d4");
    d5_param     = base->GetParamInt("d5");
    d6_param     = base->GetParamInt("d6");
    d7_param     = base->GetParamInt("d7");
    pnmodes_param= base->GetParamInt("pnmodes");
    ppid_param   = base->GetParamInt("ppid");

    nalias_param = base->GetParamInt("nalias");
    lalias_param = base->GetParamInt("lalias");
    

    //now loop over particles and normalize BRs
    //This is very important for the mass sampling
    Int_t pids=0;
    Int_t key=makeDataBase()->GetEntry("std_set");
    if (key<0) Warning("PStaticData","std_set not found");
    Int_t listkey=-1;
    while (makeDataBase()->MakeListIterator(key, "snpart","slink" , &listkey)) {
	//loop over all particles
	NormParticleBRbyKey(listkey);
	pids++;
    }

    //This is for the ctor of the PChannelModel
    defkey_param = base->MakeParamInt("defkey", "Model Def Key");

    //Add some useful aliases
    //no idea if this is the correct place, it the moment it should be OK
    AddAlias("dilepton","dielectron");
    AddAlias("g","gamma");

    AddAlias("w","omega");
    AddAlias("D0","Delta0");
    AddAlias("D+","Delta+");
    AddAlias("D++","Delta++");
    AddAlias("D-","Delta-");
    AddAlias("D0","Delta(1232)0");
    AddAlias("D+","Delta(1232)+");
    AddAlias("D++","Delta(1232)++");
    AddAlias("D-","Delta(1232)-");


    AddAlias("DP330","Delta(1600)0");
    AddAlias("DP33+","Delta(1600)+");
    AddAlias("DP33++","Delta(1600)++");
    AddAlias("DP33-","Delta(1600)-");

    AddAlias("DS310","Delta(1620)0");
    AddAlias("DS31+","Delta(1620)+");
    AddAlias("DS31++","Delta(1620)++");
    AddAlias("DS31-","Delta(1620)-");

    AddAlias("NP11+","N*(1440)+");
    AddAlias("NP110","N*(1440)0");
    AddAlias("ND13+","N*(1520)+");
    AddAlias("ND130","N*(1520)0");
    AddAlias("NS11+","N*(1535)+");
    AddAlias("NS110","N*(1535)0");

    //Ingo's additions:
    makeDataBase()->SetParamDouble ("D+", "lmass", 1.); 
    //--> lmass=0 leads to endless loop in PHadronDecayM1:: maxBWWeight
    makeDataBase()->SetParamDouble ("D+", "umass", 3.);

    system_alloc_verbosity = GetBatchValue("_system_alloc_verbosity");
    (*system_alloc_verbosity)=1.;

    *GetBatchValue("_system_weight_version") = 1.;
    *GetBatchValue("_system_unstable_width") = 0.0001;

    *GetBatchValue("_event_plane")           = 0.;
    *GetBatchValue("_event_impact_param")    = 0.;
    *GetBatchValue("_event_vertex_x")        = 0.;
    *GetBatchValue("_event_vertex_y")        = 0.;
    *GetBatchValue("_event_vertex_z")        = 0.;
    *GetBatchValue("_event_seqnr")           = -1.;


    Double_t *version= GetBatchValue("_system_version");

    if (!strncmp("CVS",version_string,3)) {
	(*version)=999.;
    } else {
	sscanf(version_string,"%lf",version);
    }
}


int PStaticData::MakeDirectoryEntry(const char * name, const char * n, const char * l, const char *ename) {
    // Adds an entry of the form "dir"
    // if entry/dir is existing, the key w/o warning is returned
    // n and l are the data base columns of the list

    Int_t skey=-1;

    if (makeDataBase()->GetParamInt(n) < 0) makeDataBase()->MakeParamInt(n,"number of links");
    if (makeDataBase()->GetParamInt(l) < 0) makeDataBase()->MakeParamInt(l,"links");

    skey=makeDataBase()->GetEntry(name);
    if (skey < 0) skey=makeDataBase()->AddEntry(name);
    if (skey < 0) return -1; //failed

    skey=makeDataBase()->GetEntry(ename);
    if (skey < 0) skey=makeDataBase()->AddEntry(ename);
    if (skey < 0) return -1; //failed
    
    Int_t *dummy;
    if (!makeDataBase()->GetParamInt(skey,l,&dummy)) {
	//fresh entry
	skey = makeDataBase()->AddListEntry(name,n,l,ename);
    }    
    return skey;

}


Double_t *PStaticData::GetBatchValue(const char * name, Int_t make_val) {
    Int_t key_a = MakeDirectoryEntry("batch_objects",NBATCH_NAME,LBATCH_NAME,name);
    //Check if double is existing
    Double_t *val;
    Int_t batch_value_param = makeDataBase()->GetParamDouble("batch_value");
    if (batch_value_param<0) 
	batch_value_param = makeDataBase()->MakeParamDouble("batch_value", "Double storage for batch");
    //cout << key_a << ":" << makeDataBase()->GetParamDouble(key_a ,batch_value_param,&val) << endl;
    if (!makeDataBase()->GetParamDouble(key_a ,batch_value_param,&val)) {
	if (make_val) {
	    Double_t * delme =  new Double_t(0.);
	    makeDataBase()->SetParamDouble(key_a,"batch_value", delme);
	    return delme;
	} else return NULL;
    }
    return val;
};



int PStaticData::GetParticleID(const char * id, int warn) {
    // pid by name
    if (!id) return 0;
    if (! makeDataBase()->GetParamInt (makeDataBase()->GetEntry((char*)id), pid_param,&i_result)) {
	Int_t key = GetAliasParent(id);
	if (key<0) {
	    if (warn) {
		Warning("GetParticleID","%s not found",id);
	    }
	    return 0;
	}
	return GetParticleIDByKey(key);
    }
    return *i_result;
};

const char * PStaticData::GetParticleName(const int & id) {
    // name by pid
    Int_t key=-1;
    if ((key=makeDataBase()->GetEntryInt(pid_param, id)) <0) {
	Warning("GetParticleName","id %i not found",id);
	return "";
    }
    makeDataBase()->GetParamString (key, name_param, &c_result);
    return c_result;
};

int PStaticData::GetParticleIDByKey(int key) {
    // pid by key
    // -1 if unvalid
    if (! makeDataBase()->GetParamInt (key, pid_param,&i_result)) 
	return -1;
    return *i_result;
}

int PStaticData::GetParticleKey(const int & id) {
    // data base key by pid
    Int_t key=-1;
    if ((key=makeDataBase()->GetEntryInt(pid_param, id)) <0) {
	Warning("GetParticleKey","id %i not found",id);
    }
    return key;
}

int PStaticData::GetDecayKey(const int & id) {
    // data base key by pid
    Int_t key=-1;
    if ((key=makeDataBase()->GetEntryInt(didx_param, id)) <0) {
	Warning("GetDecayKey","id %i not found",id);
    }
    return key;
}

int PStaticData::IsParticle(const int & id, const char * name) {
    // does pid correspond to given name?
    if (strcmp(name,GetParticleName(id))==0) return 1;
    return 0;
};

int PStaticData::IsParticleValid(const int & id) { 
    // check id range by id
    // Returns a "1" in any case if data base is filled with id
    // PFireball particles return 0;

    Int_t key=-1;
    if ((key=makeDataBase()->GetEntryInt(pid_param, id)) >=0) {
	return id; //if in data base in any case say yes
    }
 
    if ((id>500) && (id<1000)) return id;

    return 0;
}
 
int PStaticData::IsParticleValid(const char * n) {
    // check id range by name

    if (!n) return 0;
    int pid = GetParticleID(n,0);
    if (pid<0) return 0;

    return IsParticleValid(pid);
}

int PStaticData::AddAlias(const char * old_name ,const char *new_name) { 
//adds an alias to primary_key, ret value is alias key
    PDataBase * base=makeDataBase();
    Int_t pkey=-1;
    
    if ((pkey=base->AddListEntry(old_name,"nalias","lalias",new_name))<0) {
	Warning("AddAlias","Name %s not found",old_name);
	return -1;
    }
    
    return pkey;
}

int PStaticData::GetAliasParent(const char *alias_name) {
    Int_t listkey=-1,*dummy;
    Int_t key = makeDataBase()->GetEntry(alias_name);
    if (key<0) return -1;
    while (makeDataBase()->MakeListIterator(key,-1,lalias_param , &listkey)) {
	if (makeDataBase()->GetParamInt (listkey, nalias_param,&dummy)) {
	    //count found
	    return listkey;
	}
    }
    return -1;
}

int PStaticData::GetAliasParent(int key) {
    Int_t listkey=-1,*dummy;
    if (key<0) return -1;
    if (!makeDataBase()->GetParamInt (key, lalias_param,&dummy)) return -1;
    while (makeDataBase()->MakeListIterator(key,-1,lalias_param , &listkey)) {
	if (makeDataBase()->GetParamInt (listkey, nalias_param,&dummy)) {
	    //count found
	    return listkey;
	}
    }
    return -1;
}
    
int PStaticData::GetSecondaryKey(int key, int defkey) {
    Int_t listkey=-1;

    Int_t *dummy;
    if (!makeDataBase()->GetParamInt (key, nalias_param,&dummy)) {
	//Warning("MakeListIterator","count %s not found",count);
	//try to find parent
	//cout << "orig key=" << key;
	key = GetAliasParent(key);
	//cout << " now: " << key << endl;
	if (key<0)
	    return -1; //avoid messages 
    }

    while (makeDataBase()->MakeListIterator(key, nalias_param,lalias_param , &listkey)) {

	//cout << "listkey " << listkey << endl;

	//loop over all secondaries
	if (makeDataBase()->GetParamInt (listkey , defkey_param, &dummy)) {
	    if (*dummy  == defkey)
		return listkey;
	}
    }
    return -1;
}


int PStaticData::AddParticle(int pid, const char * name, double mass) {
    clearFreezeOut();
    PDataBase * base=makeDataBase();
    Int_t pkey=-1;
    
    //check if pid already exists
    if (pid>=0)
	if (base->GetEntryInt("pid",pid) >=0 ) {
	    Warning("AddParticle","pid %i already used in data base",pid);
	    return -1;
	}

    if (base->GetEntry(name) >=0) {
	Warning("AddParticle","Name %s already used in data base",name);
	return -1;
    }

    if (pid<0) {
	pid=0;
	while (makeDataBase()->GetEntryInt("pid",pid) >=0 ) {
	    pid++;
	}
	if (*system_alloc_verbosity)
	    Info("AddParticle","(%s) PID: %i", PRINT_AUTO_ALLOC,pid);

    }

    if ((pkey=base->AddListEntry("std_set","snpart","slink",name))<0) {
	Warning("AddParticle","particle header not found");
	return -1;
    }

    Int_t *ii=new int(pid);  //never destructed, but called only once!
    if (!base->SetParamInt (pkey, "pid", ii)) {
	delete ii;
	return -1;
    }
    Double_t  *dd=new double(0.);
    if (!base->SetParamDouble (pkey, "width", dd)) {
	delete dd;
	return -1;
    }
    dd=new Double_t(mass);
    if (!base->SetParamDouble (pkey, "mass", dd)) { 
	delete dd;
	return -1;
    }
    dd=new Double_t(mass);
    if (!base->SetParamDouble (pkey, "ethreshold", dd)) { 
	delete dd;
	return -1;
    }
    ii=new int(0);
    if (!base->SetParamInt (pkey, "meson", ii))
	return -1;
    ii=new int(0);
    if (!base->SetParamInt (pkey, "baryon", ii))
	return -1;
    ii=new int(0); //Meson by default: No limitation in LMass
    if (!base->SetParamInt (pkey, "lepton", ii))
	return -1;
    ii=new int(0);
    if (!base->SetParamInt (pkey, "charge", ii))
	return -1;
    ii=new int(0);
    if (!base->SetParamInt (pkey, "spin", ii))
	return -1;
    ii=new int(0);
    if (!base->SetParamInt (pkey, "parity", ii))
	return -1;
    ii=new int(0);
    if (!base->SetParamInt (pkey, "ispin", ii))
	return -1;
    ii=new int(0);
    if (!base->SetParamInt (pkey, "pythiakf", ii))
	return -1;
    ii=new int(0);  
    if (!base->SetParamInt (pkey, "widx", ii))
	return -1;
    ii=new int(0);
    if (!base->SetParamInt (pkey, "tdepth", ii))
	return -1;
    ii=new int(0);
    if (!base->SetParamInt (pkey, "hdepth", ii))
	return -1;
    return pid;
}


void  PStaticData::PrintParticle(int pid) {
    return PrintParticleByKey(makeDataBase()->GetEntryInt("pid",pid));
};


void  PStaticData::PrintParticleByKey(int key) {
    makeDataBase()->ListEntries(key,0,"name,pid,width,mass");

    Int_t listkey=-1, alias_printed=0;
    Int_t *dummy;
    if (makeDataBase()->GetParamInt (key,nalias_param,&dummy)) {
	while (makeDataBase()->MakeListIterator(key, "nalias", "lalias", &listkey)) {
	    alias_printed=1;
	    cout << "Alias: " << makeDataBase()->GetName(listkey) << endl;;
	}
    }
    if (alias_printed) cout << endl;

    if (!GetParticleNChannelsByKey(key)) {
	cout << "This particle is stable" << endl;
	return;
    }

    cout << "This particle decays via the following modes:" << endl;

    //now loop over decay modes
    listkey=-1;

    while (makeDataBase()->MakeListIterator(key, "pnmodes", "link", &listkey)) {
	PrintDecayByKey(listkey);
    }
  
};

int PStaticData::GetParticleKF(const int id) {
    // return Pythia6 kf code
    if (!id) return 0;
    if (! makeDataBase()->GetParamInt (pid_param, id, pkf_param, &i_result)) {
	Warning("GetParticleKF","id %i not found",id);
	return 0;
    }
    return *i_result;
}

int PStaticData::GetParticleIDbyKF(const int kf) {
    // return Id corresponding to Pythia6 kf code
    if (!kf) return 0;
    int *id,key;
    
    if ((key=makeDataBase()->GetEntryInt(pkf_param, kf)) <0) {
	Warning("GetParticleIDbyKF","invalid kf code %i",*id);
	return 0;
    }
    if (!makeDataBase()->GetParamInt(key,pid_param,&id))
	return 0;
    return *id;
}


int PStaticData::IsParticleMeson(const int & id) {
    // is meson?, by pid
    if (!id) return 0;
    if (! makeDataBase()->GetParamInt (pid_param, id, meson_param, &i_result)) {
	// if not found obviously not a meson
	return 0;
    }
    return *i_result;
};

int PStaticData::IsParticleHadron(const int & id) {
    // is hadron?, by pid
    if (IsParticleMeson(id) || GetParticleBaryon(id)) return 1;
    return 0;
};

int PStaticData::GetParticleBaryon(const int & id) {
    // baryon number by pid
    if (! makeDataBase()->GetParamInt (pid_param, id, baryon_param, &i_result)) {
	// if not found obviously not a baryon
	return 0;
    }
    return *i_result;
};

void PStaticData::SetParticleBaryon(const char * id, Int_t num) {
    // baryon number by name
    if (! makeDataBase()->GetParamInt (pid_param,  GetParticleID((char*)id), baryon_param, &i_result)) {
	// if not found, set param
	makeDataBase()->SetParamInt ((char*)id, "baryon", num);
    }
    *i_result=num;
}

void  PStaticData::SetParticleMeson(const char * id, Int_t num) {
 // set meson number
    if (! makeDataBase()->GetParamInt (pid_param,  GetParticleID((char*)id), meson_param, &i_result)) {
	// if not found, set param
	makeDataBase()->SetParamInt ((char*)id, "meson", num);
    }
    *i_result=num;
}

void  PStaticData::SetParticleLepton(const char * id, Int_t num) {
    if (! makeDataBase()->GetParamInt (pid_param,  GetParticleID((char*)id), lepton_param, &i_result)) {
	// if not found, set param
	makeDataBase()->SetParamInt ((char*)id, "lepton", num);
    }
    *i_result=num;
}

void  PStaticData::SetParticleCharge(const char * id, Int_t charge) {
    if (! makeDataBase()->GetParamInt (pid_param,  GetParticleID((char*)id), charge_param, &i_result)) {
	// if not found, set param
	makeDataBase()->SetParamInt ((char*)id, "charge", charge);
    }
    *i_result=charge;
}

void  PStaticData::SetParticleSpin(const char * id, Int_t spin) {
    if (! makeDataBase()->GetParamInt (pid_param,  GetParticleID((char*)id), spin_param, &i_result)) {
	// if not found, set param
	makeDataBase()->SetParamInt ((char*)id, "spin", spin);
    }
    *i_result=spin;
}

void  PStaticData::SetParticleIsospin(const char * id, Int_t isospin) {
    if (! makeDataBase()->GetParamInt (pid_param,  GetParticleID((char*)id), ispin_param, &i_result)) {
	// if not found, set param
	makeDataBase()->SetParamInt ((char*)id, "ispin", isospin);
    }
    *i_result=isospin;
}

void  PStaticData::SetParticleParity(const char * id, Int_t parity) {
    if (! makeDataBase()->GetParamInt (pid_param,  GetParticleID((char*)id), parity_param, &i_result)) {
	// if not found, set param
	makeDataBase()->SetParamInt ((char*)id, "parity", parity);
    }
    *i_result=parity;
}

int PStaticData::GetParticleLepton(const int & id) {
    // lepton number by pid
    if (! makeDataBase()->GetParamInt (pid_param, id, lepton_param, &i_result)) {
	// if not found obviously not a lepton
	return 0;
    }
    return *i_result;
};

int PStaticData::GetParticleCharge(const int & id) {
    // charge by pid
    if (! makeDataBase()->GetParamInt (pid_param, id, charge_param, &i_result)) {
	cout << "PStaticData::Charge: " << id << " not found" << endl;
	return 0;
    }
    return *i_result;
};

int PStaticData::GetParticleCharge(const char * id) {
    // charge by name
    if (! makeDataBase()->GetParamInt (pid_param, GetParticleID(id), charge_param ,&i_result)) {
	cout << "PStaticData::Charge: " << id << " not found" << endl;
	return 0;
    }
    return *i_result;
};

int PStaticData::GetParticleSpin(const int & id) {
    // 2 x J by pid
    if (! makeDataBase()->GetParamInt (pid_param, id, spin_param, &i_result)) {
	cout << "PStaticData::GetParticleSpin: " << id << " not found" << endl;
	return 0;
    }
    return *i_result;
};

int PStaticData::GetParticleSpin(const char * id) {
    // 2 x J by name
    if (! makeDataBase()->GetParamInt ( pid_param, GetParticleID(id),spin_param ,&i_result)) {
	cout << "PStaticData::GetParticleSpin: " << id << " not found" << endl;
	return 0;
    }
    return *i_result;
};

int PStaticData::GetParticleIsospin(const int & id) {
    // 2 x I by pid
    if (! makeDataBase()->GetParamInt (pid_param, id, ispin_param, &i_result)) {
	cout << "PStaticData::GetParticleIsospin: " << id << " not found" << endl;
	return 0;
    }
    return *i_result;
};

int PStaticData::GetParticleIsospin(const char * id) {
    // 2 x I by name
    if (! makeDataBase()->GetParamInt (pid_param, GetParticleID(id),ispin_param , &i_result)) {
	cout << "PStaticData::GetParticleIsospin: " << id << " not found" << endl;
	return 0;
    }
    return *i_result;
};

int PStaticData::GetParticleParity(const int & id) {
    // parity (0 if irrelevant)
    if (! makeDataBase()->GetParamInt (pid_param, id, parity_param, &i_result)) {
	cout << "PStaticData::GetParticleParity: " << id << " not found" << endl;
	return 0;
    }
    return *i_result;
    
};

int PStaticData::GetParticleParity(const char * id) {
    // parity (0 if irrelevant)
    if (! makeDataBase()->GetParamInt (pid_param, GetParticleID(id), parity_param, &i_result)) {
	cout << "PStaticData::GetParticleParity: " << id << " not found" << endl;
	return 0;
    }
    return *i_result;
};

double PStaticData::GetParticleMass(const int & id) {
    // mass by id

    if (! makeDataBase()->GetParamDouble (pid_param, id, mass_param, &d_result)) {
	cout << "PStaticData::GetParticleMass: " << id << " not found" << endl;
	return 0;
    }
    return *d_result;
};

double PStaticData::GetParticleMass(const char * id) {
    // mass by name
    if (! makeDataBase()->GetParamDouble (pid_param, GetParticleID(id),mass_param, &d_result)) {
	cout << "PStaticData::GetParticleMass: " << id << " not found" << endl;
	return 0;
    }
    return *d_result;
};

double PStaticData::GetParticleMassByKey(const int & id) {
    // mass by key
    if (! makeDataBase()->GetParamDouble (id, mass_param, &d_result)) {
	cout << "PStaticData::GetParticleMassByKey: " << id << " not found" << endl;
	return 0;
    }
    return *d_result;
};

void PStaticData::SetParticleMass(Int_t id, Float_t mass) {
    //reset mass
    if (! makeDataBase()->GetParamDouble (pid_param, id,mass_param, &d_result)) {
	cout << "PStaticData::SetParticleMass: " << id << " not found" << endl;
    }
    *d_result=mass;
}

void PStaticData::SetParticleMass(const char * id, Float_t mass) {
    //reset mass
    if (! makeDataBase()->GetParamDouble (pid_param,GetParticleID(id) ,mass_param, &d_result)) {
	cout << "PStaticData::SetParticleMass: " << id << " not found" << endl;
    }
    *d_result=mass;
}


const char *PStaticData:: GetDecayName(Int_t id) {
    if (! makeDataBase()->GetParamString (didx_param, id ,name_param, &c_result)) {
	return "NONAME";
    }
    return c_result;
}

const char *PStaticData:: GetDecayNameByKey(Int_t key) {
    if (! makeDataBase()->GetParamString (key,name_param, &c_result)) {
	return "NONAME";
    }
    return c_result;
}

void PStaticData:: SetParticleTotalWidth(const char * id, Float_t wid) {   // set new total width
    if (! makeDataBase()->GetParamDouble (pid_param, GetParticleID(id), width_param, &d_result)) {
	cout << "PStaticData::SetTWidth: " << id << " not found" << endl;
    }
    *d_result=wid;

    SetParticleEmin(GetParticleID(id),GetParticleMass(id)-2*wid);

//TODO
//   Int_t n=PNModes[id], p0=GetPosition(id), p;
//   for (p=p0; p<p0+n; ++p) if (PWidx[p]!=-1) PWidx[p]=0;
//   TWidx[id]=0;
    //printf("\n*** Warning:  Use only at start of macro! ***\n\n");
}

void PStaticData:: SetParticleTotalWidth(Int_t id, Float_t wid) {   // set new total width
    if (! makeDataBase()->GetParamDouble (pid_param, id, width_param, &d_result)) {
	cout << "PStaticData::SetTWidth: " << id << " not found" << endl;
    }
    *d_result=wid;
//TODO
//   Int_t n=PNModes[id], p0=GetPosition(id), p;
//   for (p=p0; p<p0+n; ++p) if (PWidx[p]!=-1) PWidx[p]=0;
//   TWidx[id]=0;
    printf("\n*** Warning:  Use only at start of macro! ***\n\n");
}



int PStaticData::GetTWidx(const int & id) {
    // width flag from index
    // 0: static width only
    // 1: use dynamic width
    //-1: Disabled
    if (! makeDataBase()->GetParamInt (pid_param, id , widx_param, &i_result)) {
	cout << "PStaticData::GetTWidx: " << id << " not found" << endl;
	return 0;
    }
    
    return *i_result;
}

int PStaticData::GetPWidx(const int & id) {
    // width flag from index
    // 0: static width only
    // 1: use dynamic width
    //-1: Disabled
    if (! makeDataBase()->GetParamInt (didx_param, id , widx_param, &i_result)) {
	cout << "PStaticData::GetPWidx: " << id << " not found" << endl;
	return 0;
    }
    return *i_result;
}

void PStaticData::SetTWidx(const int & id, const int & v) {
    // width flag from index
    // 0: static width only
    // 1: use dynamic width
    
    if (! makeDataBase()->GetParamInt (pid_param, id , widx_param, &i_result)) {
	cout << "PStaticData::SetTWidx: " << id << " not found" << endl;
    }
    *i_result=v;
}

void PStaticData::SetPWidx(const int & id, const int & v) {
    // width flag from index
    // 0: static width only
    // 1: use dynamic width
    
    if (! makeDataBase()->GetParamInt (didx_param, id , widx_param, &i_result)) {
	cout << "PStaticData::GetPWidx: " << id << " not found" << endl;
    }
    *i_result=v;
}

void PStaticData::SetTWidthMesh(const int & id, PMesh * mesh) {

    Int_t key=makeDataBase()->GetEntryInt("pid",id);
    if (! makeDataBase()->SetParamTObj (key ,"mesh", mesh)) {
	cout << "PStaticData::SetTWidx: " << id << " not possible" << endl;
    }
};


PMesh * PStaticData::GetTWidthMesh(const int & id) {
    if (! makeDataBase()->GetParamTObj (pid_param, id , mesh_param, &t_result)) {
	cout << "PStaticData::GetTWidthMesh: " << id << " not found" << endl;
	return NULL;
    }
    return (PMesh *) t_result;
};

void PStaticData::SetPWidthMesh(const int & id, PMesh * mesh) {
    Int_t key=makeDataBase()->GetEntryInt("didx",id);
    if (! makeDataBase()->SetParamTObj (key ,"mesh", mesh)) {
	cout << "PStaticData::SetPWidthMesh: " << id << " not possible" << endl;
    }
};


PMesh * PStaticData::GetPWidthMesh(const int & id) {
    if (! makeDataBase()->GetParamTObj (didx_param, id , mesh_param, &t_result)) {
	cout << "PStaticData::GetPWidthMesh: " << id << " not found" << endl;
	return NULL;
    }
    return (PMesh *) t_result;
};




void PStaticData::SetTF1(const int & id, TF1 * mesh) {
    Int_t key=makeDataBase()->GetEntryInt("pid",id);
    if (! makeDataBase()->SetParamTObj (key ,"tf1", mesh)) {
	cout << "PStaticData::SetTF1: " << id << " not possible" << endl;
    }
};


TF1 * PStaticData::GetTF1(const int & id) {
    if (! makeDataBase()->GetParamTObj (pid_param, id , tf1_param, &t_result)) {
//	cout << "PStaticData::GetTF1: " << id << " not found" << endl;
	return NULL;
    }
    return (TF1 *) t_result;
};



Double_t PStaticData:: GetDecayEmin(const int & id) {
    if (! makeDataBase()->GetParamDouble (didx_param, id , ethreshold_param, &d_result)) {
	cout << "PStaticData::GetDecayEmin: " << id << " not found" << endl;
	return 0;
    }
    return *d_result;   
}

int PStaticData:: GetDecayBRFlag(int didx){
    if (! makeDataBase()->GetParamInt (didx_param, didx, brflag_param, &i_result)) {
	return 0;
    }
    return *i_result;  
}

void PStaticData:: SetDecayBRFlag(int didx, int flag) {
    //Set the BR flag for the "didx", thus the total normalization
    //is used (when flag=1) instead of the pole mass
    if (! makeDataBase()->GetParamInt (didx_param, didx , brflag_param, &i_result)) {
	Int_t *ii=new Int_t(flag);
	if (!makeDataBase()->SetParamInt (GetDecayKey(didx), "brflag", ii)) {
	    delete ii;
	    return;
	}
    } else {
	*i_result=flag;
    }
}

void PStaticData:: SetTotalNormalization(char * p,int flag) {
    //Wrapper function which enables the total normalization
    //for a particle "p"
    Int_t listkey=-1;
    Int_t key = GetParticleKey(p);
    if (key<0) {
	Warning("SetTotalNormalization","Invalid particle");
	return;
    }
    
    while (makeDataBase()->MakeListIterator(key, "pnmodes", "link", &listkey)) {
	if (! makeDataBase()->GetParamInt (listkey, brflag_param, &i_result)) {
	    Int_t *ii=new Int_t(flag);
	    if (!makeDataBase()->SetParamInt (listkey, "brflag", ii)) {
		delete ii;
		return;
	    }
	} else {
	    *i_result=flag;
	}
    }
}

void PStaticData:: SetDecayEmin(const int & id, const double v) {
    if (! makeDataBase()->GetParamDouble (didx_param, id , ethreshold_param, &d_result)) {
	cout << "PStaticData::SetDecayEmin: " << id << " not found" << endl;
    }
    *d_result=v;   
}

Double_t PStaticData:: GetParticleEmin(const int & id) {
    if (! makeDataBase()->GetParamDouble (pid_param, id , ethreshold_param, &d_result)) {
	cout << "PStaticData::GetParticleEmin: " << id << " not found" << endl;
	return 0;
    }
    return *d_result;   
}

void PStaticData:: SetParticleEmin(const int & id, const double v) {
    if (! makeDataBase()->GetParamDouble (pid_param, id , ethreshold_param, &d_result)) {
	cout << "PStaticData::SetParticleEmin: " << id << " not found" << endl;
    }
    *d_result=v;   
}

Double_t PStaticData::GetParticleLMass(const int & id) {
    if (! makeDataBase()->GetParamDouble (pid_param, id , lmass_param, &d_result)) {
	Warning("GetParticleLMass","id %i not found",id);
	return 0;
    }
    return *d_result;   
}

void PStaticData::SetParticleLMass(const char * id, const double v) {
    if (! makeDataBase()->GetParamDouble (pid_param, GetParticleID((char*)id) , lmass_param, &d_result)) {
	// if not found, set param
	makeDataBase()->SetParamDouble ((char *)id, "lmass", v);
	return;
    }
    *d_result=v;
}

Double_t PStaticData::GetParticleUMass(const int & id) {
    if (! makeDataBase()->GetParamDouble (pid_param, id , umass_param, &d_result)) {
	Warning("GetParticleUMass","id %i not found",id);
	return 0;
    }
    return *d_result;   
}

void PStaticData::SetParticleUMass(const char * id, const double v) {
    if (! makeDataBase()->GetParamDouble (pid_param, GetParticleID((char*)id) , umass_param, &d_result)) {
	// if not found, set param
	makeDataBase()->SetParamDouble ((char *)id, "umass", v);
	return;
    }
    *d_result=v;
}

int PStaticData:: GetTDepth(const int & id) {
    if (! makeDataBase()->GetParamInt (pid_param, id , tdepth_param, &i_result)) {
	cout << "PStaticData::GetTDepth: " << id << " not found" << endl;
	return 0;
    }
    return *i_result;   
}

void PStaticData:: SetTDepth(const int & id, const int & depth) {
    if (! makeDataBase()->GetParamInt (pid_param, id , tdepth_param, &i_result)) {
	cout << "PStaticData::GetTDepth: " << id << " not found" << endl;
    }
    *i_result=depth;   
}

int PStaticData:: GetHDepth(const int & id) {
    if (! makeDataBase()->GetParamInt (pid_param, id , hdepth_param, &i_result)) {
	cout << "PStaticData::GetHDepth: " << id << " not found" << endl;
	return 0;
    }
    return *i_result;   
}

void PStaticData:: SetHDepth(const int & id, const int & depth) {
    if (! makeDataBase()->GetParamInt (pid_param, id , hdepth_param, &i_result)) {
	cout << "PStaticData::GetHDepth: " << id << " not found" << endl;
    }
    *i_result=depth;   
}

Double_t PStaticData:: GetDecayBR(Int_t id) {
    if (! makeDataBase()->GetParamDouble (didx_param, id ,br_param, &d_result)) {
	cout << "PStaticData::GetDecayBR: " << id << " not found" << endl;
    }
    return *d_result;

};

Double_t PStaticData:: GetDecayPartialWidth(Int_t id) {
    double res = GetDecayBR(id)*GetParticleTotalWidth(GetDecayParent(id));
    return res;
};

Double_t PStaticData:: GetDecayPartialWidthByKey(Int_t id) {
    double res = GetDecayBR(id)*GetParticleTotalWidth(GetDecayParentByKey(id));
    return res;
};


bool PStaticData:: SetDecayBR(int didx, double br, int mode) {
    //Mode: see SetDecayBR(const char *parent, const char * daughters)

    clearFreezeOut();
    Int_t key = GetDecayKey(didx);
    Double_t *brorig;
    if (!makeDataBase()->GetParamDouble (key, "brorig",&brorig)) return kFALSE;

    if (mode) {
	Double_t other = 1.-GetDecayBR(didx);
	br = (other * br)/(1.-br);
    }

    *brorig = br;
    
    return NormParticleBRbyKey(GetParticleKey(GetDecayParentByKey(key)));
}

bool PStaticData:: SetDecayBRByKey(int key, double br, int mode) {
    //Mode: see SetDecayBR(const char *parent, const char * daughters)
    clearFreezeOut();
    Double_t *brorig;
    if (!makeDataBase()->GetParamDouble (key, "brorig",&brorig)) return kFALSE;

    if (mode) {
	Double_t other = 1.-GetDecayBR(GetDecayIdxByKey(key));
	br = (other * br)/(1.-br);
    }

    *brorig = br;
    return NormParticleBRbyKey(GetParticleKey(GetDecayParentByKey(key)));
}

bool PStaticData:: SetDecayBR(const char *parent, const char * daughters, 
			      double br, int mode) {
    // Resets the decay branching ratio of an existing decay
    // Handle this functon with care!
    // Do not use this function to change the weighting etc.
    // because the change of the decay b.r. will affect the
    // shape of the parent resonance
    // 
    // If you have really physics arguments to change the b.r.
    // (e.g. update of PDG values)
    // you *must* use this function at the very beginning of your macro,
    // but at least before the "GetWidth" is called
    //
    // Mode flag:
    // 0: Add the new b.r. to the existing ones + re-weighting
    // 1: No re-weighting (in this case br must be <1.)


    //get parent key
    Int_t parent_key = makeDataBase()->GetEntry(parent);
    if (parent_key < 0) {
	cout << "PStaticData::SetDecayBR: parent " << parent << " not found in data base" << endl;
	return kFALSE;
    }

    Int_t pid=GetParticleIDByKey(parent_key);
    if (pid<0) {
	return kFALSE;
    }

    //now check the daughters 
    //TODO: All this should go into a help function
    char *arr1[7];
    Int_t arr1_s=7; //max 7 decay products
    
    PUtils::Tokenize(daughters, ",", arr1, &arr1_s);

    Int_t pids[8];
    pids[0] = pid;
    for (int pat = 0; pat < arr1_s; pat++) {
	pids[pat+1] = GetParticleID(arr1[pat]);
	//cout << pids[pat+1]<<endl;
    }

    Int_t decay_key = GetDecayKey(pids, arr1_s);

    if (decay_key<0) {
	Warning("SetDecayBR","Decay %s -> %s not defined",parent,daughters);
	return kFALSE;
    }
    
    return SetDecayBRByKey(decay_key, br, mode);
}

// This was the old Kagarlis method:
//   if (brn<1 || brn>PNModes[id]) return;
//   Int_t pos = GetPosition(id);
//   if (pos!=-1) {
//      Int_t idx = pos+brn-1;
//      PBR[idx] = br;
//      if (freeze) {
//        PWidx[idx] = -1;   // Set BR flag to static
//        printf("Branch %3i set static to %f\n",idx,br);
//      }
//      printf("\n*** Warning:  Use only at start of macro! ***\n\n");
//   }
//}

bool PStaticData:: NormParticleBR(Int_t id) { 
// normalize branching ratios for particle id
    //TODO: check if this works
    Int_t key=makeDataBase()->GetEntryInt(pid_param,id);
    return NormParticleBRbyKey(key);
}

bool PStaticData:: NormParticleBRbyKey(Int_t key) { 
// normalize branching ratios for particle id

    Int_t listkey=-1;
    double sum=0.;
    Double_t *br,*brorig,*width;
    Double_t twidth=GetParticleTotalWidthByKey(key);
    
    if (GetParticleNChannelsByKey(key)==0) return kFALSE;

    while (makeDataBase()->MakeListIterator(key, "pnmodes", "link", &listkey)) {
	if (!makeDataBase()->GetParamDouble (listkey, "brorig",&brorig)) return kFALSE;
	sum+=*brorig;
    }
    
    listkey=-1;
    if (sum>0.) 
	while (makeDataBase()->MakeListIterator(key, "pnmodes", "link", &listkey)) {
	    if (!makeDataBase()->GetParamDouble (listkey, "brorig",&brorig)) return kFALSE;
	    if (!makeDataBase()->GetParamDouble (listkey, "br",&br)) return kFALSE;
	    *br = *brorig/sum;
	    if (makeDataBase()->GetParamDouble (listkey, "width",&width)) {
		*width= *br * twidth;
	    }
	}
    return kTRUE;
}


void PStaticData::FreezeDecayBR(Int_t id, Int_t brn) { // set BR static
    Fatal("FreezeDecayBR","not implemented");

    //TODO
//   Int_t pos = GetPosition(id);
//   if (pos!=-1) {
//      Int_t idx = pos+brn-1;
//      PWidx[idx] = -1;
//      printf("Branch %3i set static\n",idx);
//      printf("\n*** Warning:  Use only at start of macro! ***\n\n");
//   }
}

int PStaticData::GetParticleNChannels(const int & id) {
    // number of decay channels by pid
    if (! makeDataBase()->GetParamInt (pid_param, id, pnmodes_param, &i_result)) {
	// if not found obviously no decay
	return 0;
    }
    return *i_result;
};

int PStaticData::GetParticleNChannels(const char * id) {
// number of decay channels by name
    if (! makeDataBase()->GetParamInt ((char*)id, "pnmodes", &i_result)) {
	return 0;
    }
    return *i_result;
};

int PStaticData::GetParticleNChannelsByKey(int key) {
// number of decay channels by key
    if (! makeDataBase()->GetParamInt (key, "pnmodes", &i_result)) {
	return 0;
    }
    return *i_result;
};

int PStaticData::IsDecayHadronic(Int_t didx) {
    Int_t tid[11];
    tid[0]=10; 

    GetDecayMode(didx, tid); 
    for (int i=1;i<tid[0];i++) {
//	cout << didx << ":" << IsParticleHadron(tid[1]) << endl;
	if (!IsParticleHadron(tid[1])) return 0;
    }
    return 1;
}

int PStaticData::AddDecay(int * ipid, int n) {
    //For internal use, e.g. PChannel

    TString *decay_string= new TString("");
    TString *decay_string2=new TString(makeStaticData()->GetParticleName(ipid[0]));
    decay_string2->Append(" --> ");	
    //first I have to check that my decay_string is big enough
    for (int i=1;i<=n;i++) {
	if (ipid[i]>1000) decay_string2->Append("(");
	decay_string->Append(makeStaticData()->GetParticleName(ipid[i]));
	decay_string2->Append(makeStaticData()->GetParticleName(ipid[i]));
	if (ipid[i]>1000) decay_string2->Append(")");
	if (i!=n) decay_string->Append(" , ");	
	if (i!=n) decay_string2->Append(" + ");
	
    }
    if (makeStaticData()->AddDecay(-1, (const char*) decay_string2->Data(), 
				   makeStaticData()->GetParticleName(ipid[0]), 
				   (const char*) decay_string->Data() , 1.0 )) {
	
	//Info("AddDecay","(%s) Decay of added: %s", PRINT_AUTO_ALLOC,decay_string->Data());
    } else {
	Warning("AddDecay","Add Decay failed");
    }
    
    return makeStaticData()->GetDecayKey(ipid, n);

}

int PStaticData::AddDecay(int didx, const char * name, const char *parent, 
			  const char * daughters , double br) {
    //make new decay with decay index "didx"
    //the didx should be a free number (or set to -1 for auto-alloc)
    //parent particle is "parent", which must be existing
    //"name" is the new unique string identifier of the decay
    //daughers is an particle array of the format "p1,p2,p3,..."
    //br is the branching ratio
    //all branching ratios are re-normalized
    clearFreezeOut();
    PDataBase *base=makeDataBase();
  
    //get parent key
    Int_t parent_key = makeDataBase()->GetEntry(parent);
    if (parent_key < 0) {
	cout << "PStaticData::addDecay: parent " << parent << " not found in data base" << endl;
	return -1;
    }

    Int_t pid=GetParticleIDByKey(parent_key);
    if (pid<0) {
	return -1;
    }
    if (didx<0) {
	if (pid<1000) didx=pid*32;
	else didx=pid;
	while (makeDataBase()->GetEntryInt("didx",didx) >=0 ) {
	    didx++;
	}
	if (*system_alloc_verbosity)
	    Info("AddDecay","(%s) Decay index %i: %s", PRINT_AUTO_ALLOC,didx,name);
    }
  


    //check if didx already exists
    if (makeDataBase()->GetEntryInt("didx",didx) >=0 ) {
	cout << "PStaticData::addDecay: didx " << didx << " already used in data base" << endl;
	return -1;
    }
    Int_t dkey = makeDataBase()->AddListEntry(parent,"pnmodes", "link",name);
    if (dkey < 0) {
	Warning("AddDecay","Unable to add decay entry to data base");
	return -1;
    }

    //now add the daughters
    char *arr1[7];
    Int_t arr1_s=7; //max 7 decay products
    

    PUtils::Tokenize(daughters, ",", arr1, &arr1_s);
    //does decay already exist?

    Int_t pids[8];
    pids[0] = pid;
    // cout << pid << endl;
    for (int pat = 0; pat < arr1_s; pat++) {
	pids[pat+1] = GetParticleID(arr1[pat]);
	//cout << pids[pat+1]<<endl;
    }
    if (GetDecayKey(pids, arr1_s)>=0) {
	Warning("AddDecay","Decay already defined");
	return -1;
    }
    

    for (int pat = 0; pat < arr1_s; pat++) {
	char * partc=arr1[pat];
	//	Int_t *pkey = new int(base->GetEntry(partc));
	//This takes into account aliases:
	Int_t *pkey = new int( makeStaticData()->GetParticleKey(partc));
	if (*pkey<0) {
	    cout << "Error: processing decay: do not find pid " << partc << endl;
	} 
	const char *ds;
	if (pat==0) ds="d1";
	if (pat==1) ds="d2";
	if (pat==2) ds="d3";
	if (pat==3) ds="d4";
	if (pat==4) ds="d5";
	if (pat==5) ds="d6";
	if (pat==6) ds="d7";
	if (pat>6) {
	    Warning("AddDecay","Too many particles");
	    return -1;
	}
	
	base->SetParamInt (dkey,ds,pkey);
    }
    //set default values
    Int_t *ii;
    makeDataBase()->GetParamInt(parent_key,"pid",&ii);
    base->SetParamInt (dkey, "ppid",ii); //set parent id
    
    ii=new int(didx);
    base->SetParamInt (dkey, "didx",ii); //decay mode index
    
    ii=new int(0);  //never destructed, but called only once!
    if (!base->SetParamInt (dkey, "widx", ii))
	return -1;
    
    Double_t * dd=new double(br);
    if (!base->SetParamDouble (dkey, "br", dd))
	return -1;
    if (!base->SetParamDouble (dkey, "brorig", new Double_t(br)))
	return -1;
    dd=new Double_t(0.);
    if (!base->SetParamDouble (dkey, "ethreshold", dd))
	return -1;
    dd=new Double_t(1.);
    if (!base->SetParamDouble (dkey, "scfactor", dd))
	return -1;

    if (pid>1000) //Do not norm compostite
	return didx;

    if (! NormParticleBRbyKey(parent_key) ) return -1;
    return didx;
    
}


void  PStaticData::PrintDecayByKey(int key) {
    Int_t products=GetDecayNProductsByKey(key);
    //  cout << products << endl;
    if (products==1)
	makeDataBase()->ListEntries(key,0,"name,didx,br,d1:name");
    if (products==2)
	makeDataBase()->ListEntries(key,0,"name,didx,br,d1:name,d2:name");
    if (products==3)
	makeDataBase()->ListEntries(key,0,"name,didx,br,d1:name,d2:name,d3:name");
    if (products==4)
	makeDataBase()->ListEntries(key,0,"name,didx,br,d1:name,d2:name,d3:name,d4:name");
    if (products==5)
	makeDataBase()->ListEntries(key,0,"name,didx,br,d1:name,d2:name,d3:name,d4:name,d5:name");
    if (products==6)
	makeDataBase()->ListEntries(key,0,"name,didx,br,d1:name,d2:name,d3:name,d4:name,d5:name,d6:name");
    if (products==7)
	makeDataBase()->ListEntries(key,0,"name,didx,br,d1:name,d2:name,d3:name,d4:name,d5:name,d6:name,d7:name");
};

int PStaticData::GetDecayNProductsByKey(const int & key) {
// retrieve number of products by mode key
    if (makeDataBase()->GetParamInt (key, d7_param, &i_result))
	return 7;
    if (makeDataBase()->GetParamInt (key, d6_param, &i_result))
	return 6;
    if (makeDataBase()->GetParamInt (key, d5_param, &i_result))
	return 5;
    if (makeDataBase()->GetParamInt (key, d4_param, &i_result))
	return 4;
    if (makeDataBase()->GetParamInt (key, d3_param, &i_result))
	return 3;
    if (makeDataBase()->GetParamInt (key, d2_param, &i_result))
	return 2;
    if (makeDataBase()->GetParamInt (key, d1_param, &i_result))
	return 1;
    return 0;
};

int PStaticData::GetDecayNProducts(const int & id) {
// retrieve number of products by mode index 
    if (makeDataBase()->GetParamInt (didx_param, id, d7_param, &i_result))
	return 7;
    if (makeDataBase()->GetParamInt (didx_param, id, d6_param, &i_result))
	return 6;
    if (makeDataBase()->GetParamInt (didx_param, id, d5_param, &i_result))
	return 5;
    if (makeDataBase()->GetParamInt (didx_param, id, d4_param, &i_result))
	return 4;
    if (makeDataBase()->GetParamInt (didx_param, id, d3_param, &i_result))
	return 3;
    if (makeDataBase()->GetParamInt (didx_param, id, d2_param, &i_result))
	return 2;
    if (makeDataBase()->GetParamInt (didx_param, id, d1_param, &i_result))
	return 1;
    return 0;
};

int PStaticData::GetDecayNProducts(const char * id){
// number of products by name     
    if (makeDataBase()->GetParamInt ((char*)id, "d7", &i_result))
	return 7;
    if (makeDataBase()->GetParamInt ((char*)id, "d6", &i_result))
	return 6;
    if (makeDataBase()->GetParamInt ((char*)id, "d5", &i_result))
	return 5;
    if (makeDataBase()->GetParamInt ((char*)id, "d4", &i_result))
	return 4;
    if (makeDataBase()->GetParamInt ((char*)id, "d3", &i_result))
	return 3;
    if (makeDataBase()->GetParamInt ((char*)id, "d2", &i_result))
	return 2;
    if (makeDataBase()->GetParamInt ((char*)id, "d1", &i_result))
	return 1;
    return 0;

}

int PStaticData::GetDecayParent(const int & id) {
    if (! makeDataBase()->GetParamInt (didx_param, id, ppid_param, &i_result)) {
	// if not found obviously no parent
	return 0;
    }
    return *i_result;
}

int PStaticData::GetDecayParentByKey(const int & id) {
    if (! makeDataBase()->GetParamInt (id, ppid_param, &i_result)) {
	// if not found obviously no parent
	return 0;
    }
    return *i_result;
}

int PStaticData:: GetDecayIdxByKey(int key) {
    // return value is -1 on failure

    Int_t * pos;
    if (!makeDataBase()->GetParamInt (key, didx_param , &pos))
	return -1;
    return *pos;
}

int PStaticData:: GetDecayIdx(int * pid, int n) {
    // decay-mode index from parent and product ids;
    // n is the size of the array
    // arguments: pointer to pid array of parent & products,
    //            number of products

    int i, nid[n];
    if (n<2||!pid) return -1;
    for (i=0;i<=n;++i) {
	if (!PStaticData::IsParticleValid(pid[i])) {
	    cout << "PStaticData:: GetDecayIdx: id not found:" << pid[i] << endl;
	    return -2;
	}
	else if (i) nid[i-1]=pid[i];        // make own copy of product array
    }
    PUtils::isort(nid,n);                 // sort the array entries

    int id=pid[0], *nm;        // parent id, number of decay modes
    //Get count info
  
    if (! makeDataBase()->GetParamInt (pid_param, id, count_param, &nm)) {
	return -3;
    }
    if (!nm) return -4;
    Int_t count;
    //now loop over decay modes
    Int_t key=makeDataBase()->GetEntryInt(pid_param, id);
    Int_t listkey=-1;
    Int_t tid[11];
    while (makeDataBase()->MakeListIterator(key, "pnmodes", "link", &listkey)) {
	tid[0]=10; 
	GetDecayModeByKey(listkey,tid);       // retrieve current mode info
	if (tid[0]==n) {                      // number of products matches input
	    PUtils::isort(tid+1,n);           // sort the current mode product pid array
	    count=0;                          // reset the match counter
	    for (i=0;i<n;++i) count+=(tid[i+1]==nid[i]);
	    if (count==n) {                   // input matched
		return GetDecayIdxByKey(listkey);          
	    }
	}
    }
  
    return -1;
}

int PStaticData:: GetDecayKey(int * pid, int n) {
    // decay-mode key from parent and product ids;
    // n is the size of the daughters array
    // arguments: pointer to pid array of parent & products,
    //            number of products

    int i, nid[n];
    if (n<1||!pid) return -1;
    for (i=0;i<=n;++i) {
	if (!PStaticData::IsParticleValid(pid[i])) {
	    Warning("GetDecayKey","id %i not found", pid[i]);
	    return -2;
	} else if (i) { 
	    nid[i-1]=pid[i];        // make own copy of product array
	}
    }

    PUtils::isort(nid,n);                 // sort the array entries

    int id=pid[0], *nm;        // parent id, number of decay modes
    //Get count info
    if (! makeDataBase()->GetParamInt (pid_param, id, count_param, &nm)) {
	return -3;
    }
    if (!nm) return -4;
    Int_t count;
    //now loop over decay modes
    Int_t key=makeDataBase()->GetEntryInt(pid_param, id);
    Int_t listkey=-1;
    Int_t tid[11];
    while (makeDataBase()->MakeListIterator(key, "pnmodes", "link", &listkey)) {
	tid[0]=10; 
	GetDecayModeByKey(listkey,tid);       // retrieve current mode info
	if (tid[0]==n) {                      // number of products matches input
	    PUtils::isort(tid+1,n);           // sort the current mode product pid array
	    count=0;                          // reset the match counter
	    for (i=0;i<n;++i) count+=(tid[i+1]==nid[i]);
	    if (count==n) {                   // input matched
		return listkey;          
	    }
	}
    }
  
    return -1;
}


void PStaticData:: GetDecayMode(int idx, int * id) {
    // Retrieves the decay mode info for a given channel n.
    // _______________________________________________________________________
    // ARGUMENTS: 1. decay-mode index (PPosition),
    // 2. pointer to existing integer array (must have large 
    // enough dimension), 3. pointer to existing double (optional).
    // _______________________________________________________________________
    // RETURNS: int array with N, i_1, i_2,..., i_N, where N is
    // the number of products for the decay mode with index n,
    // followed by the pids of these products
    //
    // N should contain the maximum number of the array-size of before
    // calling this function
    // if the maximum is exceeded, or any other error occurs, N is set to 0
    Int_t *d1_key,*d2_key,*d3_key,*d4_key,*d5_key, 
	*d6_key,*d7_key,found=0;;
    if (makeDataBase()->GetParamInt (didx_param, idx, d1_param, &d1_key)) 
	found=1;
    if (makeDataBase()->GetParamInt (didx_param, idx, d2_param, &d2_key)) 
	found=2;
    if (makeDataBase()->GetParamInt (didx_param, idx, d3_param, &d3_key)) 
	found=3;
    if (makeDataBase()->GetParamInt (didx_param, idx, d4_param, &d4_key)) 
	found=4;
    if (makeDataBase()->GetParamInt (didx_param, idx, d5_param, &d5_key)) 
	found=5;
    if (makeDataBase()->GetParamInt (didx_param, idx, d6_param, &d6_key)) 
	found=6;
    if (makeDataBase()->GetParamInt (didx_param, idx, d7_param, &d7_key)) 
	found=7;
 
    if (found>*id) {
	cout << "PStaticData:: size too low" << endl;
	*id=0;
	return;
    } else *id=found;
    Int_t *d1_p=NULL,*d2_p=NULL,*d3_p=NULL,*d4_p=NULL,*d5_p=NULL,*d6_p=NULL,*d7_p=NULL;
    if ((found>0)&& !makeDataBase()->GetParamInt (*d1_key,pid_param,&d1_p)) {
	cout << "PStaticData:: Mode: unable to unpack key1 " << *d1_key << endl;
	*id=0;	
    }
    else if ((found>1)&& !makeDataBase()->GetParamInt (*d2_key,pid_param,&d2_p)) {
	cout << "PStaticData:: Mode: unable to unpack key2 " << *d2_key << endl;
	*id=0;
    } 
    else if ((found>2)&& !makeDataBase()->GetParamInt (*d3_key,pid_param,&d3_p)) {
	cout << "PStaticData:: Mode: unable to unpack key3 " << *d3_key << endl;
	*id=0;
    } 
    else if ((found>3)&& !makeDataBase()->GetParamInt (*d4_key,pid_param,&d4_p)) {
	cout << "PStaticData:: Mode: unable to unpack key4 " << *d4_key << endl;
	*id=0;
    } 
    else if ((found>4)&& !makeDataBase()->GetParamInt (*d5_key,pid_param,&d5_p)) {
	cout << "PStaticData:: Mode: unable to unpack key5 " << *d5_key << endl;
	*id=0;
    } 
    else if ((found>5)&& !makeDataBase()->GetParamInt (*d6_key,pid_param,&d6_p)) {
	cout << "PStaticData:: Mode: unable to unpack key6 " << *d6_key << endl;
	*id=0;
    } 
    else if ((found>6)&& !makeDataBase()->GetParamInt (*d7_key,pid_param,&d7_p)) {
	cout << "PStaticData:: Mode: unable to unpack key7 " << *d7_key << endl;
	*id=0;
    } 
    if (d1_p) id[1]=*d1_p;
    if (d2_p) id[2]=*d2_p;
    if (d3_p) id[3]=*d3_p;
    if (d4_p) id[4]=*d4_p;
    if (d5_p) id[5]=*d5_p;
    if (d6_p) id[6]=*d6_p;
    if (d7_p) id[7]=*d7_p;
}


void PStaticData:: GetDecayModeByKey(int idx, int * id) {
    //same as above, but by data base key    
    Int_t *d1_key,*d2_key,*d3_key,*d4_key,*d5_key,
	*d6_key,*d7_key,found=0;
    if (makeDataBase()->GetParamInt (idx, d1_param, &d1_key)) 
	found=1;
    if (makeDataBase()->GetParamInt (idx, d2_param, &d2_key)) 
	found=2;
    if (makeDataBase()->GetParamInt (idx, d3_param, &d3_key)) 
	found=3;
    if (makeDataBase()->GetParamInt (idx, d4_param, &d4_key)) 
	found=4;
    if (makeDataBase()->GetParamInt (idx, d5_param, &d5_key)) 
	found=5;
    if (makeDataBase()->GetParamInt (idx, d6_param, &d6_key)) 
	found=6;
    if (makeDataBase()->GetParamInt (idx, d7_param, &d7_key)) 
	found=7;

    if (found>*id) {
	cout << "PStaticData:: size too low" << endl;
	*id=0;
	return;
    } else *id=found;
    //   cout << "id is " << *id << endl;
    Int_t *d1_p=NULL,*d2_p=NULL,*d3_p=NULL,*d4_p=NULL,*d5_p=NULL,*d6_p=NULL,*d7_p=NULL;
    if (*id>0 && !makeDataBase()->GetParamInt (*d1_key,pid_param,&d1_p)) {
	cout << "PStaticData:: ModeByKey: unable to unpack key1 " << *d1_key << endl;
	*id=0;
    }
    else if (*id>1 && !makeDataBase()->GetParamInt (*d2_key,pid_param,&d2_p)) {	
	cout << "PStaticData:: ModeByKey: unable to unpack key2 " << *d2_key << endl;
	*id=0;
    } 
    else if (*id>2 && !makeDataBase()->GetParamInt (*d3_key,pid_param,&d3_p)) {
	cout << "PStaticData:: ModeByKey: unable to unpack key3 " << *d3_key << endl;
	*id=0;
    }
    else if (*id>3 && !makeDataBase()->GetParamInt (*d4_key,pid_param,&d4_p)) {
	cout << "PStaticData:: ModeByKey: unable to unpack key4 " << *d4_key << endl;
	*id=0;
    } 
    else if (*id>4 && !makeDataBase()->GetParamInt (*d5_key,pid_param,&d5_p)) {
	cout << "PStaticData:: ModeByKey: unable to unpack key5 " << *d5_key << endl;
	*id=0;
    }
    else if (*id>5 && !makeDataBase()->GetParamInt (*d6_key,pid_param,&d6_p)) {
	cout << "PStaticData:: ModeByKey: unable to unpack key6 " << *d6_key << endl;
	*id=0;
    } 
    else if (*id>6 && !makeDataBase()->GetParamInt (*d7_key,pid_param,&d7_p)) {
	cout << "PStaticData:: ModeByKey: unable to unpack key7 " << *d7_key << endl;
	*id=0;
    }
    if (d1_p) id[1]=*d1_p;
    if (d2_p) id[2]=*d2_p;
    if (d3_p) id[3]=*d3_p;
    if (d4_p) id[4]=*d4_p;
    if (d5_p) id[5]=*d5_p;
    if (d6_p) id[6]=*d6_p;
    if (d7_p) id[7]=*d7_p;

}


double PStaticData::GetParticleTotalWidth(const int & id) {
    // PWidth[id]
    if (! makeDataBase()->GetParamDouble (pid_param, id, width_param, &d_result)) {
	cout << "PStaticData::GetParticleTotalWidth: " << id << " not found" << endl;
	return 0;
    }
    return *d_result;
}

double PStaticData::GetParticleTotalWidthByKey(const int & key) {
    // PWidth[id]
    if (! makeDataBase()->GetParamDouble (key, width_param, &d_result)) {
	cout << "PStaticData::GetParticleTotalWidthByKey: " << key << " not found" << endl;
	return 0;
    }
    return *d_result;
}


void listParticle(int id = -1) {
   
    if (id<0) {
        makeStaticData(); //this fill data base
	makeDataBase()->ListEntries(-1,1,"pid,*name");
    }
    else
        makeStaticData()->PrintParticle(id);
}

void listModes(int id = -1) {
   makeStaticData()->PrintParticle(id);
}
 


ClassImp(PStaticData)
