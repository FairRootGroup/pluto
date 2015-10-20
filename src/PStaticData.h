// Author: M.A. Kagarlis
// Written: 31.01.99
// Revised: (IF)
// PDataUtil Class Header

#ifndef _PSTATICDATA_H_
#define _PSTATICDATA_H_

#include "TROOT.h"
#include "TF1.h"
#include <iostream>

#include "PDataBase.h"
#include "PMesh.h"


#define MAX_DAUGHTERS 7

using namespace std;

class PStaticData;
PStaticData* makeStaticData();
PStaticData& fStaticData();


void listParticle(int);
void listParticle(const char * id);
void listModes(int id);
void listModes(const char * id);



class PStaticData : public TObject {
      
 public:

    PStaticData();

    //Freeze out -> Important to check this
    void SetFreezeOut(void) {
	freeze=kTRUE;
    };

    void clearFreezeOut(void) {
	//re-loop PStdModels
	freeze=kFALSE;
    };
    
    bool GetFreezeOut(void) {
	return freeze;
    }

    //Common stuff
    int   AddAlias(const char * old_name ,
		   const char *new_name); //adds an alias, ret value is alias key
    int   GetAliasParent(const char *alias_name);
    int   GetAliasParent(int key);
    
    int   MakeDirectoryEntry(const char * name, const char * n, const char * l,
			     const char * ename);  //adds a directory, ret value is key
    Double_t *GetBatchValue(const char * name, Int_t make_val=1);
    int   GetSecondaryKey(int key, int defkey); //Loops over the alias entries and look for matching defkey

    //Particle get methods
    int   GetParticleID(const char * id, int warn=1);     // pid by name
    int   GetParticleIDByKey(int key);        // pid by key
    const char *GetParticleName(const int & id);  // name by pid
    int   GetParticleKey(const int & id);     // data base key by pid
    int   GetParticleKey(const char * id) {
	return GetParticleKey(GetParticleID(id));
    };
    int   IsParticle(const int & id, const char * name); // does pid correspond to given name?
    int   IsParticleValid(const int & id);          // check id range by id
    int   IsParticleValid(const char * n);    // check id range by name

    int   AddParticle(int pid, const char * name, double mass);
    void  PrintParticle(int pid);
    void  PrintParticleByKey(int pid);
    void  PrintParticle(const char * id) {
	PrintParticle(GetParticleID(id));
    };

    int   GetParticleKF(const int Id);     // return Pythia6 kf code
    int   GetParticleIDbyKF(const int kf); // return Id corresponding to Pythia6 kf code

    int   IsParticleMeson(const int & id);   // is meson?, by pid
    void  SetParticleMeson(const char* id, Int_t num=1);   // set meson number

    int   IsParticleHadron(const int & id);  // is hadron?, by pid
    int   GetParticleBaryon(const int & id);    // baryon number by pid
    void  SetParticleBaryon(const char * id, Int_t num=1);    // baryon number by name

    int   GetParticleLepton(const int & id);    // lepton number by pid
    void  SetParticleLepton(const char * id, Int_t num=1);

    int   GetParticleCharge(const int & id);    // charge by pid
    int   GetParticleCharge(const char * id);   // charge by name
    void  SetParticleCharge(const char * id, Int_t charge);

    int   GetParticleSpin(const int & id);      // 2 x J by pid
    int   GetParticleSpin(const char * id);     // 2 x J by name
    void  SetParticleSpin(const char * id, Int_t spin);

    int   GetParticleIsospin(const int & id);   // 2 x I by pid
    int   GetParticleIsospin(const char * id);  // 2 x I by name
    void  SetParticleIsospin(const char * id, Int_t isospin);

    int   GetParticleParity(const int & id);    // parity (0 if irrelevant)
    int   GetParticleParity(const char * id);   // parity (0 if irrelevant)
    void  SetParticleParity(const char * id, Int_t parity);

    double GetParticleMass(const int & id);   // mass by id  
    double GetParticleMass(const char * id);  // mass by name
    double GetParticleMassByKey(const int & id);
    void SetParticleMass(Int_t id, Float_t mass); //reset mass
    void SetParticleMass(const char * id, Float_t mass); //reset mass

    int GetParticleNChannels(const int & id); // number of decay channels by pid  
    int GetParticleNChannels(const char * id);// number of decay channels by name
    int GetParticleNChannelsByKey(int id);

    void   SetParticleTotalWidth(Int_t id, Float_t wid);
    void   SetParticleTotalWidth(const char * id, Float_t wid);
    double GetParticleTotalWidth(const int & id); // -->PWidth[id]
    double GetParticleTotalWidth(const char * id) {
	return GetParticleTotalWidth(GetParticleID(id));
	};
    double GetParticleTotalWidthByKey(const int & id); 

    double GetParticleEmin(const int & id); // Returns the energy threshold for the particle
    void   SetParticleEmin(const int & id, const double v);

    //Change particle range:
    double GetParticleLMass(const int & id); // Returns lower mass (used in all samplings)
    double GetParticleLMass(const char * id) {
	return GetParticleLMass(GetParticleID(id));
    };
    double GetParticleUMass(const int & id); // Returns upper mass (used in all samplings)
    double GetParticleUMass(const char * id) {
	return GetParticleUMass(GetParticleID(id));
    };
    void SetParticleLMass(const int & id, const double v) { // Set lower mass (used in all samplings)
	SetParticleLMass(GetParticleName(id),v);
    };
    void SetParticleLMass(const char * id, const double v); // Set lower mass (used in all samplings)

    void SetParticleUMass(const int & id, const double v) { // Set upper mass (used in all samplings)
	SetParticleUMass(GetParticleName(id),v);
    };
    void SetParticleUMass(const char * id, const double v);// Set upper mass (used in all samplings)


 
    bool NormParticleBR(Int_t id); // normalize branching ratios for particle id
    bool NormParticleBRbyKey(Int_t key);
    void SetTotalNormalization(char * p,int flag=1);

    //Decay methods
    void FreezeDecayBR(Int_t id, Int_t brn); // Set BR  (BUGBUG->brn nomally unknown)
    bool SetDecayBR(int didx, double br, int mode);
    bool SetDecayBR(const char *parent, const char * daughters , double br, int mode);
    bool SetDecayBRByKey(int key, double br, int mode);

    Double_t GetDecayBR(Int_t id);
    Double_t GetDecayPartialWidth(Int_t id);
    Double_t GetDecayPartialWidthByKey(Int_t id);
    const char *GetDecayName(Int_t id);
    const char *GetDecayNameByKey(Int_t key);

    Int_t IsDecayHadronic(Int_t didx);
    

    int  AddDecay(int didx, const char * name, const char *parent, 
		  const char * daughters , double br);
    int  AddDecay(const char * name, const char *parent, 
		  const char * daughters , double br) {
	return AddDecay(-1, name, parent, 
			daughters , br); };
    int  AddDecay(int * ipid, int n);
    void PrintDecayByKey(int key);

    int GetDecayNProducts(const int &); // retrieve number of products by mode index 
    int GetDecayNProducts(const char *);// number of products by name 
    int GetDecayNProductsByKey(const int & key);

    int GetDecayParent(const int &); // parent pid from decay mode index
    int GetDecayParentByKey(const int &); 
    void GetDecayMode(const int, int * n); // retrieve product number, pids, by mode index
    void GetDecayModeByKey(const int, int * n); 

    int GetDecayIdx(int * pid, int n); // decay-mode index from parent and product ids; ->getChannel
    int GetDecayKey(int * pid, int n); 
    int GetDecayKey(const int & id);
    int GetDecayIdxByKey(int key);
     
    int GetDecayBRFlag(int didx);
    void SetDecayBRFlag(int didx, int flag);
     

    double GetDecayEmin(const int & idx); // Returns the energy threshold for the decay mode with index=idx
    void   SetDecayEmin(const int & idx, const double v);

    //for dynamic stuff
    int GetTWidx(const int &); // total width flag from index
    int GetPWidx(const int &); // partial width flag from index
    void SetTWidx(const int &, const int &); // total width flag from index
    void SetPWidx(const int &, const int &); // partial width flag from index

    int GetTDepth(const int &); 
    void SetTDepth(const int &, const int &); 
    int GetHDepth(const int &); 
    void SetHDepth(const int &, const int &); 

    void SetTWidthMesh(const int &, PMesh * mesh);
    PMesh * GetTWidthMesh(const int &);

    void SetPWidthMesh(const int &, PMesh * mesh);
    PMesh * GetPWidthMesh(const int &);

    void SetTF1(const int &, TF1 * mesh);
    TF1 * GetTF1(const int &);

     
    // global friend functions
    //Keeping for backward compatibility
  
    friend void listParticle(int id=-1);
    // list particles in data base and their properties, by particle pid

    friend void listParticle(const char * id) {
	listParticle(makeStaticData()->GetParticleID(id));
    };
    // list particles in data base and their properties, by particle code name
    
    friend void listModes(int id=-1);
    // list decay modes in data base, by particle pid
    
    friend void listModes(const char * id){
	listModes(makeStaticData()->GetParticleID(id));
    };
    // list decay modes in data base, by particle name
    
    

 private:
    Int_t *i_result;
    const char *c_result;
    Double_t   *d_result;
    TObject    *t_result;
    Int_t  pid_param;
    Int_t  name_param;
    Int_t  meson_param,baryon_param,lepton_param;
    Int_t  charge_param, spin_param, ispin_param;
    Int_t  parity_param, mass_param, width_param;
    Int_t  pkf_param, didx_param;
    Int_t  widx_param, mesh_param, tf1_param, ethreshold_param, lmass_param, umass_param;
    Int_t  tdepth_param, hdepth_param,br_param,brorig_param,count_param;
    Int_t  d1_param,d2_param,d3_param, pnmodes_param, ppid_param;
    Int_t  d4_param,d5_param,d6_param,d7_param;
    Int_t  brflag_param;
    Int_t  nalias_param,lalias_param, defkey_param;

    Bool_t freeze;

    Double_t  *system_alloc_verbosity;

    ClassDef(PStaticData,0) //Pluto Static Data Wrapper

};

#endif // _PDATAUTIL_H_







	
