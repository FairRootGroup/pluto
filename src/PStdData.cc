/////////////////////////////////////////////////////////////////////
//  Pluto Particle Data (The standard version)
//
//  This class contains a permanent data base of particles
//  and their properties, functions to load additional
//  temporary particles and decay modes, as well as functions
//  to build decay widths and branching ratios and sample masses
//  of arbitrary hadronic resonances.
//  The PID convention is consistent with GEANT, with up to
//  999 particles supported.
//
//                             Author:  M.A. Kagarlis
//                             Written: 31/01/99
//                             Revised: 24/05/2004 R.H.
//                             Revised particle tables: 09/06/2004 R.H.
//                             Bug fixes: 13/10/2004 R.H.
//                             mu+mu-: 02/02/2005 R.H.
//                             N1535: 25/01/2007 R.H.
//
// The static arrays habe been copied from the original PData class
// to a "filler" to have more flexibility (IF 28.4.2007)
//
// Ref 1: Particle Data Group, Review of Particle Properties,
//        Phys. Rev. D54 (1996) 1 (and earlier editions)
//////////////////////////////////////////////////////////////////////

// local arrays
#include "TArrayI.h"
#include "TArrayD.h"

#include "PDataBase.h"
#include "PStdData.h"
#include <math.h>

PStdData& fStdData()
 {
   static PStdData* ans = new PStdData();
   return *ans;
 }

PStdData * makeStdData()
 {
   return &fStdData();
 }

PStdData::PStdData() {
    Info("PStdData()","(%s)", PRINT_CTOR);
    
    disable=0;
}

PStdData::~PStdData() {
    cout << "Killing standard particle data table" << endl;
}


Bool_t PStdData::fillDataBase(void) {

    if (disable) return kTRUE;

    PPosition = new int[maxnumpar];
    resetPosition() ;

    PDataBase *base = makeDataBase();
    //make particle header
    base->MakeParamInt("setx", "Particle set index");
    base->MakeParamInt("snpart", "Number of particles");
    base->MakeParamInt("slink", "Particle link");

    Int_t pidkey=base->MakeParamInt("pid", "Pluto particle ID");
    
    base->SetFastKey(pidkey,100000);
    base->MakeParamDouble("width", "Particle static width [GeV]");
    base->MakeParamDouble("mass",  "Particle pole mass [GeV]");
    base->MakeParamInt("meson", "Meson number");
    base->MakeParamInt("baryon", "Baryon number");
    base->MakeParamInt("lepton", "Lepton number");
    base->MakeParamInt("charge", "Charge");
    base->MakeParamInt("spin", "Spin");
    base->MakeParamInt("parity", "Parity");
    base->MakeParamInt("ispin", "Isospin");
    base->MakeParamInt("pythiakf", "Pythia KF code");

    base->MakeParamInt("pnmodes", "Number of decay modes");
    base->MakeParamInt("link", "Decay link");
    base->MakeParamInt("d1", "Decay product 1");
    base->MakeParamInt("d2", "Decay product 2");
    base->MakeParamInt("d3", "Decay product 3");
    base->MakeParamInt("d4", "Decay product 4");
    base->MakeParamInt("d5", "Decay product 5");
    base->MakeParamInt("d6", "Decay product 6");
    base->MakeParamInt("d7", "Decay product 7");
    base->MakeParamInt("ppid", "Parent Pluto ID"); //parent pid
    Int_t idxkey=base->MakeParamInt("didx", "Decay index");
    base->MakeParamDouble("br", "Branching ratio");
    base->MakeParamDouble("brorig", "Original branching ratio (for norm.)");
    base->SetFastKey(idxkey,1000*32);

    base->MakeParamInt("widx", "Width flag"); //partial/total width switched on by default
    base->MakeParamInt("tdepth", "total depth");
    base->MakeParamInt("hdepth", "hadronic depth");
    base->MakeParamDouble("ethreshold", "Energy threshold");
    base->MakeParamDouble("scfactor",  "Self consistency factor");
    base->MakeParamInt("sccount", "Self consistency count (max tries)");
    

    base->MakeParamTObj("mesh", "Mesh object");
    base->MakeParamTObj("tf1", "TF1 object");
    base->MakeParamTObj("model", "Basic PModel");

    base->MakeParamInt("maxmesh", "Number of Mesh points");
    base->MakeParamDouble("lmass", "Lower mass");
    base->MakeParamDouble("umass", "Upper mass");

    base->MakeParamInt("brflag", "=1: use integral normalization");

    //alias
    base->MakeParamInt("nalias", "Number of aliases");
    base->MakeParamInt("lalias", "Alias link");

    //this is for the decay manager
    base->MakeParamTObj("decaychannel", "Decay channel for decay manager");

    //...for the decay_all option
    base->MakeParamTObj("stackchannel", "Stack of channels for the decayall option in PReaction");

    Int_t skey=-1;
    if ((skey=base->AddEntry("std_set"))<0) return kFALSE;


    int pkey=-1,*ii;
    Double_t *dd;
    for (int i=0; i<maxnumpar; i++) {
//	if ((pkey=base->AddEntry(PStdData::PName[i]))<0) return kFALSE;
//	Int_t dkey = base->AddListEntry(PStdData::PName[i],"pnmodes", "link",PMDescription [j]);
	if ((pkey=base->AddListEntry("std_set","snpart","slink",PStdData::PName[i]))<0) return kFALSE;
	ii=new int(i);  //never destructed, but called only once!
	if (!base->SetParamInt (pkey, "pid", ii))
	    return kFALSE;
	if (!base->SetParamDouble (pkey, "width", &(PStdData::PWidth[i])))
	    return kFALSE;
	if (!base->SetParamDouble (pkey, "mass", &(PStdData::PMass[i])))
	    return kFALSE;
	if (PMeson[i])
	    if (!base->SetParamInt (pkey, "meson", &(PStdData::PMeson[i])))
		return kFALSE;
	if (PBaryon[i])
	    if (!base->SetParamInt (pkey, "baryon", &(PStdData::PBaryon[i])))
		return kFALSE;
	if (PLepton[i])
	    if (!base->SetParamInt (pkey, "lepton", &(PStdData::PLepton[i])))
		return kFALSE;
	if (!base->SetParamInt (pkey, "charge", &(PStdData::PCharge[i])))
	    return kFALSE;
	if (!base->SetParamInt (pkey, "spin", &(PStdData::PJ[i])))
	    return kFALSE;
	if (!base->SetParamInt (pkey, "parity", &(PStdData::PParity[i])))
	    return kFALSE;
	if (!base->SetParamInt (pkey, "ispin", &(PStdData::PI[i])))
	    return kFALSE;
	if (!base->SetParamInt (pkey, "pythiakf", &(PStdData::Pkf[i])))
	    return kFALSE;
	ii=new int(0);  //never destructed, but called only once!
	//0 means on, -1 means off
	if (!base->SetParamInt (pkey, "widx", ii))
	    return kFALSE;
	ii=new int(0);
	if (!base->SetParamInt (pkey, "tdepth", ii))
	    return kFALSE;
	ii=new int(0);
	if (!base->SetParamInt (pkey, "hdepth", ii))
	    return kFALSE;
	dd=new Double_t(PStdData::PMass[i]-2*PStdData::PWidth[i]); //BUGBUG look to ethreshold later -> only 1st guess!!!
	if (!base->SetParamDouble (pkey, "ethreshold", dd))
	    return kFALSE;

	//Adding Fireballs!
	char * name = new char[100];
	sprintf(name,"Fireball: %s",PStdData::PName[i]);
	if ((pkey=base->AddListEntry("std_set","snpart","slink",name))<0) return kFALSE;
	ii=new int(i+500);  //never destructed, but called only once!
	if (!base->SetParamInt (pkey, "pid", ii))
	    return kFALSE;
	if (!base->SetParamDouble (pkey, "mass", &(PStdData::PMass[i])))
	    return kFALSE;
	if (!base->SetParamDouble (pkey, "width", &(PStdData::PWidth[i])))
	    return kFALSE;
    }
    
    for (int i=0; i<maxnumpar; i++) {
//	Int_t pkey=base->getEntry(PStdData::PName[i]);
	//now the decays
	int jmin=0, jmax=PStdData::maxnummodes;
	jmin=PPosition[i];
	jmax=jmin+PStdData::PNModes[i];
	//debug info
//	cout << "decay of " << PStdData::PName[i] << endl;

	for (int j=jmin;j<jmax;++j) {

	    TString s=PMode[j];         // retrieve decay mode string

	    Int_t dkey = base->AddListEntry(PStdData::PName[i],"pnmodes", "link",PMDescription [j]);

	    ii=new int(i);
	    base->SetParamInt (dkey, "ppid",ii); //set parent id

	    ii=new int(j+1); //Shifted by 1 because 0 is not good
	    base->SetParamInt (dkey, "didx",ii); //decay mode index

	    dd=new Double_t(0.); //BUGBUG look to ethreshold
	    if (!base->SetParamDouble (dkey, "ethreshold", dd))
		return kFALSE;

	    ii=new int(0);  //never destructed, but called only once!
	    if (!base->SetParamInt (dkey, "widx", ii))
		return kFALSE;
	    
	    if (!base->SetParamDouble (dkey, "br", &(PStdData::PBR[j])))
		return kFALSE;

	    if (!base->SetParamDouble (dkey, "brorig", new Double_t(PStdData::PBR[j])))
		return kFALSE;

	    dd=new Double_t(0.);
	    if (!base->SetParamDouble (dkey, "ethreshold", dd))
		return kFALSE;
	    
	    dd=new Double_t(0.); //partial width resetted later
	    if (!base->SetParamDouble (dkey, "width", dd))
		return kFALSE;
	    
	    dd=new Double_t(1.); //sc factor=1.
	    if (!base->SetParamDouble (dkey, "scfactor", dd))
		return kFALSE;

	    Int_t len=s.Length();
	    Int_t res=len%3;
	    Int_t np= len/3+(res>0);
	    Int_t pid, k=0;
	    for (int ii=0;ii<np;++ii) {   // loop over product particles
		Int_t m=(!ii&&res)?res:3;        // number of digits in current pid
		pid=0;              // reset id
		for (int jj=1;jj<=m;++jj) pid += (*(s(k+jj-1,1).Data())-48)
				       *int(pow(10.,m-jj));
		//cout << s << endl;
		//cout << pid << endl;
		Int_t *pkey = new int(base->GetEntryInt("pid",pid));
		if (*pkey<0) {
		    cout << "Error: processing decay: do not find pid " << pid << endl;
		}
		if (ii==0) base->SetParamInt (dkey, "d1",pkey);
		if (ii==1) base->SetParamInt (dkey, "d2",pkey);
		if (ii==2) base->SetParamInt (dkey, "d3",pkey);
		if (ii==3) base->SetParamInt (dkey, "d4",pkey);
		if (ii==4) base->SetParamInt (dkey, "d5",pkey);
		if (ii==5) base->SetParamInt (dkey, "d6",pkey);
		if (ii==6) base->SetParamInt (dkey, "d7",pkey);
		if (ii>6) cout << "PStdData::fillDataBase: More then 7 decay products not supported:" 
			       << PMDescription [j] << endl;
		k+=m;

	    } //end decay products

	} //end decay mode loop
    }
    
    


    disable=1;
    return kTRUE;

}


//Please do not touch the tables below

#define nmax 999     // maximum number of supported particles
#define mnpar 70     // number of particles stored permanently
#define mnmodes 220  // number of decay modes stored permanently


// particle naming convention for use with the PParticle constructor
const char *PStdData::NAME[mnpar] = {
  "dummy","g",
  "e+", "e-", "nu", "mu+", "mu-", "pi0", "pi+",
  "pi-", "K0L", "K+", "K-", "n", "p", "anti_p",
  "K0S", "eta", "Lambda", "Sigma+", "Sigma0",
  "Sigma-", "Xi0", "Xi-", "Omega", "anti_n",
  "anti_Lambda", "anti_Sigma-", "anti_Sigma0",
  "anti_Sigma+","anti_Xi0", "anti_Xi+",
  "anti_Omega+","File", "D0", "D++", "D+", "D-", "NP11+", 
  "ND13+", "NS11+", "rho0", "rho+", "rho-",
  "BOZO", "d", "t", "alpha", "BOZO2", "He3", //BOZO->BOZO2
  "dimuon", "dilepton", "w", "eta'", "sigma",
  "phi", "DP330", "DP33++", "DP33+","DP33-",
  "DS310", "DS31++", "DS31+","DS31-",
  "NP110", "ND130", "NS110", "J/Psi", "Psi'","pn"};

// particle masses (GeV/c**2)
//const double MASS[mnpar]={
double PStdData::MASS[mnpar]={
  /* 0: dummy     */  0.0,            /* 1: Photon    */  0.0,
  /* 2: Positron  */  0.51099906e-3,  /* 3: Electron  */  0.51099906e-3,
  /* 4: Neutrino  */  0.0,            /* 5: mu+       */  0.105658389,
  /* 6: mu-       */  0.105658389,    /* 7: pi0       */  0.1349764,
  /* 8: pi+       */  0.13956995,     /* 9: pi-       */  0.13956995,
  /*10: K0 long   */  0.497672,       /*11: K+        */  0.493677,
  /*12: K-        */  0.493677,       /*13: Neutron   */  0.93956563,
  /*14: Proton    */  0.93827231,     /*15: Antiproton*/  0.93827231,
  /*16: K0 short  */  0.497672,       /*17: Eta       */  0.54745,
  /*18: Lambda    */  1.115684,       /*19: Sigma+    */  1.18937,
  /*20: Sigma0    */  1.19255,        /*21: Sigma-    */  1.197436,
  /*22: Xi0       */  1.3149,         /*23: Xi-       */  1.32132,
  /*24: Omega     */  1.67245,        /*25: Antineutrn*/  0.93956563,
  /*26: Antilambda*/  1.115684,       /*27: Antisigma-*/  1.18937, //1.197436,
  /*28: Antisigma0*/  1.19255,        /*29: Antisigma+*/  1.197436, //1.18937,
  /*30: Antixi0   */  1.3149,         /*31: Antixi+   */  1.32132,
  /*32: Antiomega+*/  1.67245,        /*33: File      */  0.0,
  /*34: Delta0    */  1.232,          /*35: Delta++   */  1.232,
  /*36: Delta+    */  1.232,          /*37: Delta-    */  1.232,
  /*38: NP11+     */  1.44,           /*39: ND13+     */  1.520,
  /*40: NS11+     */  1.535,          /*41: rho0      */  0.7699,
  /*42: rho+      */  0.7699,         /*43: rho-      */  0.7699,
  /*44: NULL      */  0.0,            /*45: Deuteron  */  1.875613,
  /*46: Tritium   */  2.80925,        /*47: Alpha     */  3.727417,
  /*48: NULL      */  0.0,            /*49: He3       */  2.80923,
  /*50: dimuon    */  0.21131678,     /*51: dilepton  */  0.001022,
  /*52: omega     */  0.78194,        /*53: eta'      */  0.9577,
  /*54: sigma     */  0.6,            /*55: phi       */  1.019413,
  /*56: Delta0*P33*/  1.6,            /*57: Delta++ *P33*/1.6,
  /*58: Delta+*P33*/  1.6,            /*59: Delta- *P33 */1.6,
  /*60: Delta0*S31*/  1.62,           /*61: Delta++ *S31*/1.62,
  /*62: Delta+*S31*/  1.62,           /*63: Delta- *S31 */1.62,
  /*64: NP110     */  1.44,           /*65: ND130     */  1.520,
  /*66: NS110     */  1.535,          /*67: J/Psi     */  3.09688,
  /*68: Psi'      */  3.68596,        /*69: pn        */  2.65

};

// particle widths (GeV/c**2)
//const double WIDTH[mnpar]={
double PStdData::WIDTH[mnpar]={
  /* 0: dummy     */  0.0,             /* 1: Photon    */  0.0,
  /* 2: Positron  */  0.0,             /* 3: Electron  */  0.0,
  /* 4: Neutrino  */  0.0,             /* 5: mu+       */  hbar/2.19703e-6,
  /* 6: mu-       */  PStdData::hbar/2.19703e-6, /* 7: pi0       */  PStdData::hbar/8.4e-17,
  /* 8: pi+       */  PStdData::hbar/2.6033e-8,  /* 9: pi-       */  PStdData::hbar/2.6033e-8,
  /*10: K0 long   */  PStdData::hbar/5.17e-8,    /*11: K+        */  PStdData::hbar/1.2386e-8,
  /*12: K-        */  PStdData::hbar/1.23861e-8, /*13: Neutron   */  PStdData::hbar/887.,
  /*14: Proton    */  0.0,             /*15: Antiproton*/  0.0,
  /*16: K0 short  */  PStdData::hbar/8.927e-11,  /*17: eta       */  1.29e-6,
  /*18: Lambda    */  PStdData::hbar/2.632e-10,  /*19: Sigma+    */  PStdData::hbar/7.99e-11,
  /*20: Sigma0    */  PStdData::hbar/7.4e-20,    /*21: Sigma-    */  PStdData::hbar/1.479e-10,
  /*22: Xi0       */  PStdData::hbar/2.9e-10,    /*23: Xi-       */  PStdData::hbar/1.639e-10,
  /*24: Omega-    */  PStdData::hbar/8.22e-11,   /*25: Antineutrn*/  PStdData::hbar/887.,
  /*26: Antilambda*/  PStdData::hbar/2.632e-10,  /*27: Antisigma-*/  PStdData::hbar/7.99e-11,
  /*28: Antisigma0*/  PStdData::hbar/7.4e-20,    /*29: Antisigma+*/  PStdData::hbar/1.479e-10,
  /*30: Antixi0   */  PStdData::hbar/2.9e-10,    /*31: Antixi+   */  PStdData::hbar/1.639e-10,
  /*32: Antiomega+*/  PStdData::hbar/8.22e-11,   /*33: File      */  0.0,
  /*34: Delta0    */  0.12,            /*35: Delta++   */  0.12,
  /*36: Delta+    */  0.12,            /*37: Delta-    */  0.12,
  /*38: NP11+     */  0.35,            /*39: ND13+     */  0.12,
  /*40: NS11+     */  0.15,            /*41: rho0      */  0.1507,
  /*42: rho+      */  0.1507,          /*43: rho-      */  0.1507,
  /*44: NULL      */  0.0,             /*45: Deuteron  */  0.0,
  /*46: Tritium   */  0.0,             /*47: Alpha     */  0.0,
  /*48: NULL      */  0.0,             /*49: He3       */  0.0,
  /*50: dimuon    */  0.0,             /*51: dilepton  */  0.0,
  /*52: omega     */  0.00843,         /*53: eta'      */  0.000201,
  /*54: sigma     */  0.5,             /*55: phi       */  0.00443,
  /*56: Delta0*P33*/  0.35,            /*57: Delta++ *P33*/0.35,
  /*58: Delta+*P33*/  0.35,            /*59: Delta- *P33 */0.35,
  /*60: Delta0*S31*/  0.15,            /*61: Delta++ *S31*/0.15,
  /*62: Delta+*S31*/  0.15,            /*63: Delta- *S31 */0.15,
  /*64: NP110     */  0.35,            /*65: ND130      */ 0.12,
  /*66: NS110     */  0.15,            /*67: J/Psi      */ 0.000087,
  /*68: Psi'      */ 0.000277,         /*69: pn         */ 0.5
};

// Pythia6 KF code
const int PStdData::PYTHIAKF[mnpar]={
  /* 0: dummy     */   0,            /* 1: Photon    */  22,
  /* 2: Positron  */ -11,            /* 3: Electron  */  11,
  /* 4: Neutrino  */  12,            /* 5: mu+       */ -13,
  /* 6: mu-       */  13,            /* 7: pi0       */ 111,
  /* 8: pi+       */ 211,            /* 9: pi-       */-211,
  /*10: K0 long   */ 130,            /*11: K+        */ 321,
  /*12: K-        */-321,            /*13: Neutron   */2112,
  /*14: Proton    */2212,            /*15: Antiproton*/-2212,
  /*16: K0 short  */ 310,            /*17: eta       */ 221,
  /*18: Lambda    */3122,            /*19: Sigma+    */3222,
  /*20: Sigma0    */3212,            /*21: Sigma-    */3112,
  /*22: Xi0       */3322,            /*23: Xi-       */3312,
  /*24: Omega-    */3334,            /*25: Antineutrn*/-2112,
  /*26: Antilambda*/-3122,           /*27: Antisigma-*/-3112,
  /*28: Antisigma0*/-3212,           /*29: Antisigma+*/-3222,
  /*30: Antixi0   */-3322,           /*31: Antixi+   */-3312,
  /*32: Antiomega+*/-3334,           /*33: File      */  0,
  /*34: Delta0    */2114,            /*35: Delta++   */2224,
  /*36: Delta+    */2214,            /*37: Delta-    */1114,
  /*38: NP11+     */  0,             /*39: ND13+     */  0,
  /*40: NS11+     */  0,             /*41: rho0      */ 113,
  /*42: rho+      */ 213,            /*43: rho-      */-213,
  /*44: NULL      */  0,             /*45: Deuteron  */  0,
  /*46: Tritium   */  0,             /*47: Alpha     */  0,
  /*48: NULL      */  0,             /*49: He3       */  0,
  /*50: dimuon    */  0,             /*51: dilepton  */  0,
  /*52: omega     */ 223,            /*53: eta'      */ 331,
  /*54: sigma     */  0,             /*55: phi       */ 333,
  /*56: Delta0*P33*/  0,             /*57: Delta++ *P33*/0,
  /*58: Delta+*P33*/  0,             /*59: Delta- *P33 */0,
  /*60: Delta0*S31*/  0,             /*61: Delta++ *S31*/0,
  /*62: Delta+*S31*/  0,             /*63: Delta- *S31 */0,
  /*64: NP110     */  0,             /*65: ND130     */  0,
  /*66: NS110     */  0,             /*67: J/Psi     */ 443,
  /*68: Psi'      */  0,             /*69: pn        */  0
};

// "is Meson" flag
const int PStdData::MESON[mnpar]={
  /* 0: dummy     */  0,             /* 1: Photon    */  0,
  /* 2: Positron  */  0,             /* 3: Electron  */  0,
  /* 4: Neutrino  */  0,             /* 5: mu+       */  0,
  /* 6: mu-       */  0,             /* 7: pi0       */  1,
  /* 8: pi+       */  1,             /* 9: pi-       */  1,
  /*10: K0 long   */  1,             /*11: K+        */  1,
  /*12: K-        */  1,             /*13: Neutron   */  0,
  /*14: Proton    */  0,             /*15: Antiproton*/  0,
  /*16: K0 short  */  1,             /*17: eta       */  1,
  /*18: Lambda    */  0,             /*19: Sigma+    */  0,
  /*20: Sigma0    */  0,             /*21: Sigma-    */  0,
  /*22: Xi0       */  0,             /*23: Xi-       */  0,
  /*24: Omega-    */  0,             /*25: Antineutrn*/  0,
  /*26: Antilambda*/  0,             /*27: Antisigma-*/  0,
  /*28: Antisigma0*/  0,             /*29: Antisigma+*/  0,
  /*30: Antixi0   */  0,             /*31: Antixi+   */  0,
  /*32: Antiomega+*/  0,             /*33: File      */  0,
  /*34: Delta0    */  0,             /*35: Delta++   */  0,
  /*36: Delta+    */  0,             /*37: Delta-    */  0,
  /*38: NP11+     */  0,             /*39: ND13+     */  0,
  /*40: NS11+     */  0,             /*41: rho0      */  1,
  /*42: rho+      */  1,             /*43: rho-      */  1,
  /*44: NULL      */  0,             /*45: Deuteron  */  0,
  /*46: Tritium   */  0,             /*47: Alpha     */  0,
  /*48: NULL      */  0,             /*49: He3       */  0,
  /*50: dimuon    */  0,             /*51: dilepton  */  0,
  /*52: omega     */  1,             /*53: eta'      */  1,
  /*54: sigma     */  1,             /*55: phi       */  1,
  /*56: Delta0*P33*/  0,             /*57: Delta++ *P33*/0,
  /*58: Delta+*P33*/  0,             /*59: Delta- *P33 */0,
  /*60: Delta0*S31*/  0,             /*61: Delta++ *S31*/0,
  /*62: Delta+*S31*/  0,             /*63: Delta- *S31 */0,
  /*64: NP110     */  0,             /*65: ND130     */  0,
  /*66: NS110     */  0,             /*67: J/Psi     */  1,
  /*68: Psi'      */  1,             /*69: pn        */  0
};

// Baryon number
const int PStdData::BARYON[mnpar]={
  /* 0: dummy     */  0,             /* 1: Photon    */  0,
  /* 2: Positron  */  0,             /* 3: Electron  */  0,
  /* 4: Neutrino  */  0,             /* 5: mu+       */  0,
  /* 6: mu-       */  0,             /* 7: pi0       */  0,
  /* 8: pi+       */  0,             /* 9: pi-       */  0,
  /*10: K0 long   */  0,             /*11: K+        */  0,
  /*12: K-        */  0,             /*13: Neutron   */  1,
  /*14: Proton    */  1,             /*15: Antiproton*/ -1,
  /*16: K0 short  */  0,             /*17: eta       */  0,
  /*18: Lambda    */  1,             /*19: Sigma+    */  1,
  /*20: Sigma0    */  1,             /*21: Sigma-    */  1,
  /*22: Xi0       */  1,             /*23: Xi-       */  1,
  /*24: Omega-    */  1,             /*25: Antineutrn*/ -1,
  /*26: Antilambda*/ -1,             /*27: Antisigma-*/ -1,
  /*28: Antisigma0*/ -1,             /*29: Antisigma+*/ -1,
  /*30: Antixi0   */ -1,             /*31: Antixi+   */ -1,
  /*32: Antiomega+*/ -1,             /*33: File      */  0,
  /*34: Delta0    */  1,             /*35: Delta++   */  1,
  /*36: Delta+    */  1,             /*37: Delta-    */  1,
  /*38: NP11+     */  1,             /*39: ND13+     */  1,
  /*40: NS11+     */  1,             /*41: rho0      */  0,
  /*42: rho+      */  0,             /*43: rho-      */  0,
  /*44: NULL      */  0,             /*45: Deuteron  */  2,
  /*46: Tritium   */  3,             /*47: Alpha     */  4,
  /*48: NULL      */  0,             /*49: He3       */  3,
  /*50: dimuon    */  0,             /*51: dilepton  */  0,
  /*52: omega     */  0,             /*53: eta'      */  0,
  /*54: sigma     */  0,             /*55: phi       */  0,
  /*56: Delta0*P33*/  1,             /*57: Delta++ *P33*/1,
  /*58: Delta+*P33*/  1,             /*59: Delta- *P33 */1,
  /*60: Delta0*S31*/  1,             /*61: Delta++ *S31*/1,
  /*62: Delta+*S31*/  1,             /*63: Delta- *S31 */1,
  /*64: NP110     */  1,             /*65: ND130     */  1,
  /*66: NS110     */  1,             /*67: J/Psi     */  0,
  /*68: Psi'      */  0,             /*69: pn        */  2
};

// Lepton number
const int PStdData::LEPTON[mnpar]={
  /* 0: dummy     */  0,             /* 1: Photon    */  0,
  /* 2: Positron  */ -1,             /* 3: Electron  */  1,
  /* 4: Neutrino  */  1,             /* 5: mu+       */ -1,
  /* 6: mu-       */  1,             /* 7: pi0       */  0,
  /* 8: pi+       */  0,             /* 9: pi-       */  0,
  /*10: K0 long   */  0,             /*11: K+        */  0,
  /*12: K-        */  0,             /*13: Neutron   */  0,
  /*14: Proton    */  0,             /*15: Antiproton*/  0,
  /*16: K0 short  */  0,             /*17: eta       */  0,
  /*18: Lambda    */  0,             /*19: Sigma+    */  0,
  /*20: Sigma0    */  0,             /*21: Sigma-    */  0,
  /*22: Xi0       */  0,             /*23: Xi-       */  0,
  /*24: Omega-    */  0,             /*25: Antineutrn*/  0,
  /*26: Antilambda*/  0,             /*27: Antisigma-*/  0,
  /*28: Antisigma0*/  0,             /*29: Antisigma+*/  0,
  /*30: Antixi0   */  0,             /*31: Antixi+   */  0,
  /*32: Antiomega+*/  0,             /*33: File      */  0,
  /*34: Delta0    */  0,             /*35: Delta++   */  0,
  /*36: Delta+    */  0,             /*37: Delta-    */  0,
  /*38: NP11+     */  0,             /*39: ND13+     */  0,
  /*40: NS11+     */  0,             /*41: rho0      */  0,
  /*42: rho+      */  0,             /*43: rho-      */  0,
  /*44: NULL      */  0,             /*45: Deuteron  */  0,
  /*46: Tritium   */  0,             /*47: Alpha     */  0,
  /*48: NULL      */  0,             /*49: He3       */  0,
  /*50: dimuon    */  0,             /*51: dilepton  */  0,
  /*52: omega     */  0,             /*53: eta'      */  0,
  /*54: sigma     */  0,             /*55: phi       */  0,
  /*56: Delta0*P33*/  0,             /*57: Delta++ *P33*/0,
  /*58: Delta+*P33*/  0,             /*59: Delta- *P33 */0,
  /*60: Delta0*S31*/  0,             /*61: Delta++ *S31*/0,
  /*62: Delta+*S31*/  0,             /*63: Delta- *S31 */0,
  /*64: NP110     */  0,             /*65: ND130     */  0,
  /*66: NS110     */  0,             /*67: J/Psi     */  0,
  /*68: Psi'      */  0,             /*69: pn        */  0
};

// Particle charge
const int PStdData::CHARGE[mnpar]={
  /* 0: dummy     */  0,             /* 1: Photon    */  0,
  /* 2: Positron  */  1,             /* 3: Electron  */ -1,
  /* 4: Neutrino  */  0,             /* 5: mu+       */  1,
  /* 6: mu-       */ -1,             /* 7: pi0       */  0,
  /* 8: pi+       */  1,             /* 9: pi-       */ -1,
  /*10: K0 long   */  0,             /*11: K+        */  1,
  /*12: K-        */ -1,             /*13: Neutron   */  0,
  /*14: Proton    */  1,             /*15: Antiproton*/ -1,
  /*16: K0 short  */  0,             /*17: eta       */  0,
  /*18: Lambda    */  0,             /*19: Sigma+    */  1,
  /*20: Sigma0    */  0,             /*21: Sigma-    */ -1,
  /*22: Xi0       */  0,             /*23: Xi-       */ -1,
  /*24: Omega-    */ -1,             /*25: Antineutrn*/  0,
  /*26: Antilambda*/  0,             /*27: Antisigma-*/ -1,
  /*28: Antisigma0*/  0,             /*29: Antisigma+*/  1,
  /*30: Antixi0   */  0,             /*31: Antixi+   */  1,
  /*32: Antiomega+*/  1,             /*33: File      */  0,
  /*34: Delta0    */  0,             /*35: Delta++   */  2,
  /*36: Delta+    */  1,             /*37: Delta-    */ -1,
  /*38: NP11+     */  1,             /*39: ND13+     */  1,
  /*40: NS11+     */  1,             /*41: rho0      */  0,
  /*42: rho+      */  1,             /*43: rho-      */ -1,
  /*44: NULL      */  0,             /*45: Deuteron  */  1,
  /*46: Tritium   */  1,             /*47: Alpha     */  2,
  /*48: NULL      */  0,             /*49: He3       */  2,
  /*50: dimuon    */  0,             /*51: dilepton  */  0,
  /*52: omega     */  0,             /*53: eta'      */  0,
  /*54: sigma     */  0,             /*55: phi       */  0,
  /*56: Delta0*P33*/  0,             /*57: Delta++ *P33*/2,
  /*58: Delta+*P33*/  1,             /*59: Delta- *P33 */-1,
  /*60: Delta0*S31*/  0,             /*61: Delta++ *S31*/2,
  /*62: Delta+*S31*/  1,             /*63: Delta- *S31 */-1,
  /*64: NP110     */  0,             /*65: ND130     */  0,
  /*66: NS110     */  0,             /*67: J/Psi     */  0,
  /*68: Psi'      */  0,             /*69: pn        */  1
};

// 2 x angular momentum
const int PStdData::SPIN[mnpar]={
  /* 0: dummy     */  0,             /* 1: Photon    */  2,
  /* 2: Positron  */  1,             /* 3: Electron  */  1,
  /* 4: Neutrino  */  1,             /* 5: mu+       */  1,
  /* 6: mu-       */  1,             /* 7: pi0       */  0,
  /* 8: pi+       */  0,             /* 9: pi-       */  0,
  /*10: K0 long   */  0,             /*11: K+        */  0,
  /*12: K-        */  0,             /*13: Neutron   */  1,
  /*14: Proton    */  1,             /*15: Antiproton*/  1,
  /*16: K0 short  */  0,             /*17: eta       */  0,
  /*18: Lambda    */  1,             /*19: Sigma+    */  1,
  /*20: Sigma0    */  1,             /*21: Sigma-    */  1,
  /*22: Xi0       */  1,             /*23: Xi-       */  1,
  /*24: Omega-    */  3,             /*25: Antineutrn*/  1,
  /*26: Antilambda*/  1,             /*27: Antisigma-*/  1,
  /*28: Antisigma0*/  1,             /*29: Antisigma+*/  1,
  /*30: Antixi0   */  1,             /*31: Antixi+   */  1,
  /*32: Antiomega+*/  3,             /*33: File      */  0,
  /*34: Delta0    */  3,             /*35: Delta++   */  3,
  /*36: Delta+    */  3,             /*37: Delta-    */  3,
  /*38: NP11+     */  1,             /*39: ND13+     */  3,
  /*40: NS11+     */  1,             /*41: rho0      */  2,
  /*42: rho+      */  2,             /*43: rho-      */  2,
  /*44: NULL      */  0,             /*45: Deuteron  */  2,
  /*46: Tritium   */  1,             /*47: Alpha     */  0,
  /*48: NULL      */  0,             /*49: He3       */  1,
  /*50: dimuon    */  2,             /*51: dilepton  */  2,
  /*52: omega     */  2,             /*53: eta'      */  0,
  /*54: sigma     */  0,             /*55: phi       */  2,
  /*56: Delta0*P33*/  3,             /*57: Delta++ *P33*/3,
  /*58: Delta+*P33*/  3,             /*59: Delta- *P33 */3,
  /*60: Delta0*S31*/  1,             /*61: Delta++ *S31*/1,
  /*62: Delta+*S31*/  1,             /*63: Delta- *S31 */1,
  /*64: NP110     */  1,             /*65: ND130     */  3,
  /*66: NS110     */  1,             /*67: J/Psi     */  2,
  /*68: Psi'      */  2,             /*69: pn        */  0
};

// Parity (+/-1, 0 if irrelevant)
const int PStdData::PARITY[mnpar]={
  /* 0: dummy     */  0,             /* 1: Photon    */ -1,
  /* 2: Positron  */  0,             /* 3: Electron  */  0,
  /* 4: Neutrino  */  0,             /* 5: mu+       */  0,
  /* 6: mu-       */  0,             /* 7: pi0       */ -1,
  /* 8: pi+       */ -1,             /* 9: pi-       */ -1,
  /*10: K0 long   */ -1,             /*11: K+        */ -1,
  /*12: K-        */ -1,             /*13: Neutron   */  1,
  /*14: Proton    */  1,             /*15: Antiproton*/ -1,
  /*16: K0 short  */ -1,             /*17: eta       */ -1,
  /*18: Lambda    */  1,             /*19: Sigma+    */  1,
  /*20: Sigma0    */  1,             /*21: Sigma-    */  1,
  /*22: Xi0       */  1,             /*23: Xi-       */  1,
  /*24: Omega-    */  1,             /*25: Antineutrn*/ -1,
  /*26: Antilambda*/ -1,             /*27: Antisigma-*/ -1,
  /*28: Antisigma0*/ -1,             /*29: Antisigma+*/ -1,
  /*30: Antixi0   */ -1,             /*31: Antixi+   */ -1,
  /*32: Antiomega+*/ -1,             /*33: File      */  0,
  /*34: Delta0    */  1,             /*35: Delta++   */  1,
  /*36: Delta+    */  1,             /*37: Delta-    */  1,
  /*38: NP11+     */  1,             /*39: ND13+     */ -1,
  /*40: NS11+     */ -1,             /*41: rho0      */ -1,
  /*42: rho+      */ -1,             /*43: rho-      */ -1,
  /*44: NULL      */  0,             /*45: Deuteron  */  0,
  /*46: Tritium   */  0,             /*47: Alpha     */  0,
  /*48: NULL      */  0,             /*49: He3       */  0,
  /*50: dimuon    */ -1,             /*51: dilepton  */ -1,
  /*52: omega     */ -1,             /*53: eta'      */ -1,
  /*54: sigma     */  1,             /*55: phi       */ -1,
  /*56: Delta0*P33*/  1,             /*57: Delta++ *P33*/1,
  /*58: Delta+*P33*/  1,             /*59: Delta- *P33 */1,
  /*60: Delta0*S31*/ -1,             /*61: Delta++ *S31*/-1,
  /*62: Delta+*S31*/ -1,             /*63: Delta- *S31 */-1,
  /*64: NP110     */  1,             /*65: ND130     */ -1,
  /*66: NS110     */ -1,             /*67: J/Psi     */ -1,
  /*68: Psi'      */ -1,             /*69: pn        */  1
};

// 2 x isospin (also 0 if irrelevant)
const int PStdData::ISPIN[mnpar]={
  /* 0: dummy     */  0,             /* 1: Photon    */  0,
  /* 2: Positron  */  0,             /* 3: Electron  */  0,
  /* 4: Neutrino  */  0,             /* 5: mu+       */  0,
  /* 6: mu-       */  0,             /* 7: pi0       */  2,
  /* 8: pi+       */  2,             /* 9: pi-       */  2,
  /*10: K0 long   */  1,             /*11: K+        */  1,
  /*12: K-        */  1,             /*13: Neutron   */  1,
  /*14: Proton    */  1,             /*15: Antiproton*/  1,
  /*16: K0 short  */  1,             /*17: eta       */  0,
  /*18: Lambda    */  0,             /*19: Sigma+    */  2,
  /*20: Sigma0    */  2,             /*21: Sigma-    */  2,
  /*22: Xi0       */  1,             /*23: Xi-       */  1,
  /*24: Omega-    */  0,             /*25: Antineutrn*/  0,
  /*26: Antilambda*/  0,             /*27: Antisigma-*/  2,
  /*28: Antisigma0*/  2,             /*29: Antisigma+*/  2,
  /*30: Antixi0   */  1,             /*31: Antixi+   */  1,
  /*32: Antiomega+*/  0,             /*33: File      */  0,
  /*34: Delta0    */  3,             /*35: Delta++   */  3,
  /*36: Delta+    */  3,             /*37: Delta-    */  3,
  /*38: NP11+     */  1,             /*39: ND13+     */  1,
  /*40: NS11+     */  1,             /*41: rho0      */  2,
  /*42: rho+      */  2,             /*43: rho-      */  2,
  /*44: NULL      */  0,             /*45: Deuteron  */  0,
  /*46: Tritium   */  0,             /*47: Alpha     */  0,
  /*48: NULL      */  0,             /*49: He3       */  0,
  /*50: dimuon    */  0,             /*51: dilepton  */  0,
  /*52: omega     */  0,             /*53: eta'      */  0,
  /*54: sigma     */  0,             /*55: phi       */  0,
  /*56: Delta0*P33*/  3,             /*57: Delta++ *P33*/3,
  /*58: Delta+*P33*/  3,             /*59: Delta- *P33 */3,
  /*60: Delta0*S31*/  3,             /*61: Delta++ *S31*/3,
  /*62: Delta+*S31*/  3,             /*63: Delta- *S31 */3,
  /*64: NP110     */  1,             /*65: ND130     */  1,
  /*66: NS110     */  1,             /*67: J/Psi     */  0,
  /*68: Psi'      */  0,             /*69: pn        */  0
};

// Number of decay modes per particle
const int PStdData::NMODES[mnpar]={
  /* 0: dummy     */  0,             /* 1: Photon    */  0,
  /* 2: Positron  */  0,             /* 3: Electron  */  0,
  /* 4: Neutrino  */  0,             /* 5: mu+       */  1,
  /* 6: mu-       */  1,             /* 7: pi0       */  2,
  /* 8: pi+       */  1,             /* 9: pi-       */  1,
  /*10: K0 long   */  6,             /*11: K+        */  6,
  /*12: K-        */  6,             /*13: Neutron   */  1,
  /*14: Proton    */  0,             /*15: Antiproton*/  0,
  /*16: K0 short  */  3,             /*17: eta       */  8,
  /*18: Lambda    */  3,             /*19: Sigma+    */  0,
  /*20: Sigma0    */  0,             /*21: Sigma-    */  0,
  /*22: Xi0       */  0,             /*23: Xi-       */  0,
  /*24: Omega-    */  0,             /*25: Antineutrn*/  1,
  /*26: Antilambda*/  0,             /*27: Antisigma-*/  0,
  /*28: Antisigma0*/  0,             /*29: Antisigma+*/  0,
  /*30: Antixi0   */  0,             /*31: Antixi+   */  0,
  /*32: Antiomega+*/  0,             /*33: File      */  0,
  /*34: Delta0    */  4,             /*35: Delta++   */  1,
  /*36: Delta+    */  4,             /*37: Delta-    */  1,
  /*38: NP11+     */ 10,             /*39: ND13+     */ 10,
  /*40: NS11+     */ 13,             /*41: rho0      */  3,
  /*42: rho+      */  1,             /*43: rho-      */  1,
  /*44: NULL      */  0,             /*45: Deuteron  */  0,
  /*46: Tritium   */  0,             /*47: Alpha     */  0,
  /*48: NULL      */  0,             /*49: He3       */  0,
  /*50: dimuon    */  1,             /*51: dilepton  */  1,
  /*52: omega     */  7,             /*53: eta'      */  7,
  /*54: sigma     */  3,             /*55: phi       */  8,
  /*56: Delta0*P33*/ 10,             /*57: Delta++ *P33*/5,
  /*58: Delta+*P33*/ 10,             /*59: Delta- *P33 */5,
  /*60: Delta0*S31*/  8,             /*61: Delta++ *S31*/4,
  /*62: Delta+*S31*/  8,             /*63: Delta- *S31 */4,
  /*64: NP110     */ 10,             /*65: ND130     */ 10,
  /*66: NS110     */ 13,             /*67: J/Psi     */  8,
  /*68: Psi'      */  8,             /*69: pn        */  2
};

// Static branching ratios:
// different charge states are included explicitly
// with the appropriate isospin factors.
// Note: The static branching ratios for decay modes of
// the type 1) hadron -> hadron + hadron and 
// 2) hadron -> N + pi + pi are averages of PDG data, and
// are only included for completeness. They are actually
// calculated as functions of mass explicitly (see Width()
// and Width1()). Static BR's are only used for modes that
// are not entirely hadronic, or involve more than two
// decay products.
//const double BRR[mnmodes]={
double PStdData::BRR[mnmodes]={
  // mu+, 1 channel
  1.,            // id=5  mu+ --> e+ + neutrino + neutrino
  // mu-, 1 channel
  1.,            // id=6  mu- --> e- + neutrino + neutrino
  // pi0, 2 channels
  0.988,         // id=7  pi0 --> photon + photon
  0.012,         // id=7  pi0 --> photon + dilepton (Dalitz)
  // pi+, 1 channel
  1.,            // id=8  pi+ --> mu+ + neutrino
  // pi-, 1 channel
  1.,            // id=9  pi- --> mu- + neutrino
  // K0L, 6 channels
  0.211,         // id=10 K0 long --> pi0 + pi0 + pi0
  0.126,         // id=10 K0 long --> pi+ + pi- + pi0
  0.136,         // id=10 K0 long --> pi+ + mu- + neutrino
  0.136,         // id=10 K0 long --> pi- + mu+ + neutrino
  0.194,         // id=10 K0 long --> pi+ + e- + neutrino
  0.194,         // id=10 K0 long --> pi- + e+ + neutrino
  // K+, 6 channels
  0.635,         // id=11 K+ --> mu+ + neutrino
  0.212,         // id=11 K+ --> pi+ + pi0
  0.056,         // id=11 K+ --> pi+ + pi+ + pi-
  0.017,         // id=11 K+ --> pi+ + pi0 + pi0
  0.032,         // id=11 K+ --> pi0 + mu+ + neutrino
  0.048,         // id=11 K+ --> pi0 + e+ + neutrino
  // K-, 6 channels
  0.635,         // id=12 K- --> mu- + neutrino
  0.212,         // id=12 K- --> pi- + pi0
  0.056,         // id=12 K- --> pi- + pi- + pi+
  0.017,         // id=12 K- --> pi- + pi0 + pi0
  0.032,         // id=12 K- --> pi0 + mu- + neutrino
  0.048,         // id=12 K- --> pi0 + e- + neutrino
  // n, 1 channel
  1.,            // id=13 n --> p + e- + neutrino
  // K0S, 3 channels
  0.6851,        // id=16 K0 short --> pi+ + pi-
  0.3129,        // id=16 K0 short --> pi0 + pi0
  0.0018,        // id=16 K0 short --> pi+ + pi- + photon
  // eta, 8 channels
  0.394,         // id=17 eta --> photon + photon
  0.325,         // id=17 eta --> pi0 + pi0 + pi0
  0.226,         // id=17 eta --> pi+ + pi- + pi0
  0.0468,        // id=17 eta --> pi+ + pi- + photon
  0.006,         // id=17 eta --> photon + dilepton (Dalitz)
  0.00031,       // id=17 eta --> photon + dimuon
  6.0e-6,        // id=17 eta --> e+ + e-
  5.8e-6,        // id=17 eta --> mu+ + mu-
  // Lambda, 3 channels
  0.639,         // id=18 Lambda --> p + pi-
  0.358,         // id=18 Lambda --> n + pi0
  0.0018,        // id=18 Lambda --> n + photon
  // anti_n, 1 channel
  1.,            // id=25 anti_n --> anti_p + e+ + neutrino
  // D0, 4 channels
  0.99/3.,       // id=34 Delta0 --> p + pi-
  0.99*2./3.,    // id=34 Delta0 --> n + pi0
  0.0055,        // id=34 Delta0 --> n + photon
  0.0055/137.,   // id=34 Delta0 --> n + dilepton
  // D++, 1 channel
  1.,            // id=35 Delta++ --> p + pi+
  // D+, 4 channels
  0.99*2./3,     // id=36 Delta+ --> p + pi0
  0.99/3.,       // id=36 Delta+ --> n + pi+
  0.0055,        // id=36 Delta+ --> p + photon
  0.0055/137.,   // id=36 Delta+ --> p + dilepton
  // D-, 1 channel
  1.,            // id=36 Delta- --> n + pi-
  // NP11+, 10 channels
  0.65/3.,       // id=38 N*(1440)+ --> p + pi0
  0.65*2./3.,    // id=38 N*(1440)+ --> n + pi+
  0.25/2.,       // id=38 N*(1440)+ --> Delta++ + pi-
  0.25/3.,       // id=38 N*(1440)+ --> Delta+ + pi0
  0.25/6.,       // id=38 N*(1440)+ --> Delta0 + pi+
  0.024585/3.,   // id=38 N*(1440)+ --> p + rho0
  0.024585*2/3., // id=38 N*(1440)+ --> n + rho+
  0.075/3.,      // id=38 N*(1440)+ --> p + pi0 + pi0
  0.075*2./3.,   // id=38 N*(1440)+ --> p + pi+ + pi-
  0.000415,      // id=38 N*(1440)+ --> p + photon
  // ND13+, 10 channels
  0.55/3.,       // id=39 N*(1520)+ --> p + pi0
  0.55*2./3.,    // id=39 N*(1520)+ --> n + pi+
  0.20/2.,       // id=39 N*(1520)+ --> Delta++ + pi-
  0.20/3.,       // id=39 N*(1520)+ --> Delta+ + pi0
  0.20/6.,       // id=39 N*(1520)+ --> Delta0 + pi+
  0.20/3.,       // id=39 N*(1520)+ --> p + rho0
  0.20*2/3.,     // id=39 N*(1520)+ --> n + rho+
  0.0449/3.,     // id=39 N*(1520)+ --> p + pi0 + pi0
  0.0449*2./3.,  // id=39 N*(1520)+ --> p + pi+ + pi-
  0.0051,        // id=39 N*(1520)+ --> p + photon
  // NS11+, 13 channels
  0.46/3.,       // id=40 N*(1535)+ --> p + pi0
  0.46*2./3.,    // id=40 N*(1535)+ --> n + pi+
  0.3875,        // id=40 N*(1535)+ --> p + eta
  0.01/2.,       // id=40 N*(1535)+ --> Delta++ + pi-
  0.01/3.,       // id=40 N*(1535)+ --> Delta+ + pi0
  0.01/6.,       // id=40 N*(1535)+ --> Delta0 + pi+
  0.04/3.,       // id=40 N*(1535)+ --> p + rho0
  0.04*2/3.,     // id=40 N*(1535)+ --> n + rho+
  0.03/3.,       // id=40 N*(1535)+ --> p + pi0 + pi0
  0.03*2./3.,    // id=40 N*(1535)+ --> p + pi+ + pi-
  0.07/3.,       // id=40 N*(1535)+ --> N*(1440)+ + pi0
  0.07*2./3.,    // id=40 N*(1535)+ --> N*(1440)0 + pi+
  0.0025,        // id=40 N*(1535)+ --> p + photon
  // rho0, 3 channels
  0.9999091,     // id=41 rho0 --> pi+ + pi-
  4.67e-5,       // id=41 rho0 --> e+ + e-
  4.55e-5,       // id=41 rho0 --> mu+ + mu-
  // rho+, 1 channel
  1.,            // id=42 rho+ --> pi+ + pi0 
  // rho-, 1 channel
  1.,            // id=43 rho- --> pi- + pi0 
  // dimuon, 1 channel
  1.,            // id=50 dimuon --> mu+ + mu-
  // dilepton, 1 channel
  1.,            // id=51 dilepton --> e+ + e-
  // w, 7 channels
  0.888,         // id=52 omega --> pi+ + pi- + pi0
  0.085,         // id=52 omega --> pi0 + photon
  0.017,         // id=52 omega --> pi+ + pi-
  5.9e-4,        // id=52 omega --> pi0 + dilepton (Dalitz)
  9.6e-5,        // id=52 omega --> pi0 + dimuon
  7.07e-5,       // id=52 omega --> e+ + e-
  8.0e-5,        // id=52 omega --> mu+ + mu-
  // eta', 7 channels
  0.443,         // id=53 eta' --> pi- + pi+ + eta
  0.295,         // id=53 eta' --> rho0 + photon
  0.209,         // id=53 eta' --> pi0 + pi0 + eta
  0.0303,        // id=53 eta' --> omega + photon
  0.0212,        // id=53 eta' --> photon + photon
  0.00156,       // id=53 eta' --> pi0 + pi0 + pi0
  0.000104,      // id=53 eta' --> photon + dimuon
  // sigma, 3 channels
  0.99999,       // id=54 sigma --> pi+ + pi-
  0.000005,      // id=54 sigma --> e+ + e-
  0.000005,      // id=54 sigma --> mu+ + mu-
  // phi, 8 channels
  0.492,         // id=55 phi --> K+ + K-
  0.338,         // id=55 phi --> K0L + K0S
  0.155,         // id=55 phi --> pi+ + pi- + pi0
  0.01297,       // id=55 phi --> eta + photon
  0.00126,       // id=55 phi --> pi0 + photon
  0.000296,      // id=55 phi --> e+ + e-
  0.000287,      // id=55 phi --> mu+ + mu-
  0.000115,      // id=55 phi --> dilepton + eta
  // DP330, 10 channels
  .175/3.,       // id=56 Delta(1600)0 --> p + pi-
  .175*2./3.,    // id=56 Delta(1600)0 --> n + pi0
  .55*8./15.,    // id=56 Delta(1600)0 --> Delta+ + pi-
  .55/15.,       // id=56 Delta(1600)0 --> Delta0 + pi0
  .55*2./5.,     // id=56 Delta(1600)0 --> Delta- + pi+
  .225/3.,       // id=56 Delta(1600)0 --> p + rho-
  .225*2./3.,    // id=56 Delta(1600)0 --> n + rho0
  .0499/3.,      // id=56 Delta(1600)0 --> N(1440)+ + rho-
  .0499*2./3.,   // id=56 Delta(1600)0 --> N(1440)0 + rho0
  .0001,         // id=56 Delta(1600)0 --> n + photon
  // DP33++, 5 channels
  .175,          // id=57 Delta(1600)++ --> p + pi+
  .55*3./5.,     // id=57 Delta(1600)++ --> Delta++ + pi0
  .55*2./5.,     // id=57 Delta(1600)++ --> Delta+ + pi+
  .225,          // id=57 Delta(1600)++ --> p + rho+
  .05,           // id=57 Delta(1600)++ --> N(1440)+ + rho+
  // DP33+, 10 channels
  .175*2./3.,    // id=58 Delta(1600)+ --> p + pi0
  .175/3.,       // id=58 Delta(1600)+ --> n + pi+
  .55*2./5.,     // id=58 Delta(1600)+ --> Delta++ + pi-
  .55/15.,       // id=58 Delta(1600)+ --> Delta+ + pi0
  .55*8./15.,    // id=58 Delta(1600)+ --> Delta0 + pi+
  .225*2./3.,    // id=58 Delta(1600)+ --> p + rho0
  .225/3.,       // id=58 Delta(1600)+ --> n + rho+
  .0499*2./3.,   // id=58 Delta(1600)+ --> N(1440)+ + rho0
  .0499/3.,      // id=58 Delta(1600)+ --> N(1440)0 + rho+
  .0001,         // id=58 Delta(1600)+ --> p + photon
  // DP33-, 5 channels
  .175,          // id=59 Delta(1600)- --> n + pi-
  .55*2./5.,     // id=59 Delta(1600)- --> Delta0 + pi-
  .55*3./5.,     // id=59 Delta(1600)- --> Delta- + pi0
  .225,          // id=59 Delta(1600)- --> n + rho-
  .05,           // id=59 Delta(1600)- --> N(1440)0 + rho-
  // DS310, 8 channels
  .25/3.,        // id=60 Delta(1620)0 --> p + pi-
  .25*2./3.,     // id=60 Delta(1620)0 --> n + pi0
  .5897*8./15.,  // id=60 Delta(1620)0 --> Delta+ + pi-
  .5897/15.,     // id=60 Delta(1620)0 --> Delta0 + pi0
  .5897*2./5.,   // id=60 Delta(1620)0 --> Delta- + pi+
  .16/3.,        // id=60 Delta(1620)0 --> p + rho-
  .16*2./3.,     // id=60 Delta(1620)0 --> n + rho0
  .0003,         // id=60 Delta(1620)0 --> n + photon
  // DS31++, 4 channels
  .25,           // id=61 Delta(1620)++ --> p + pi+
  .5897*3./5.,   // id=61 Delta(1620)++ --> Delta++ + pi0
  .5897*2./5.,   // id=61 Delta(1620)++ --> Delta+ + pi+
  .16,           // id=61 Delta(1620)++ --> p + rho+
  // DS31+, 8 channels
  .25*2./3.,     // id=62 Delta(1620)+ --> p + pi0
  .25/3.,        // id=62 Delta(1620)+ --> n + pi+
  .5897*2./5.,   // id=62 Delta(1620)+ --> Delta++ + pi-
  .5897/15.,     // id=62 Delta(1620)+ --> Delta+ + pi0
  .5897*8./15.,  // id=62 Delta(1620)+ --> Delta0 + pi+
  .16*2./3.,     // id=62 Delta(1620)+ --> p + rho0
  .16/3.,        // id=62 Delta(1620)+ --> n + rho+
  .0003,         // id=62 Delta(1620)+ --> p + photon
  // DS31-, 4 channels
  .25,           // id=63 Delta(1620)- --> n + pi-
  .5897*2./5.,   // id=63 Delta(1620)- --> Delta0 + pi-
  .5897*3./5.,   // id=63 Delta(1620)- --> Delta- + pi0
  .16,           // id=63 Delta(1620)- --> n + rho-
  // NP110, 10 channels
  0.65/3.,       // id=64 N*(1440)0 --> p + pi-
  0.65*2./3.,    // id=64 N*(1440)0 --> n + pi0
  0.25/2.,       // id=64 N*(1440)0 --> Delta+ + pi-
  0.25/3.,       // id=64 N*(1440)0 --> Delta0 + pi0
  0.25/6.,       // id=64 N*(1440)0 --> Delta- + pi+
  0.024585/3.,   // id=64 N*(1440)0 --> p + rho-
  0.024585*2/3., // id=64 N*(1440)0 --> n + rho0
  0.075*2./3.,   // id=64 N*(1440)0 --> n + pi+ + pi-
  0.075/3.,      // id=64 N*(1440)0 --> n + pi0 + pi0
  0.000415,      // id=64 N*(1440)0 --> n + photon
  // ND130, 10 channels
  0.55*2./3.,    // id=65 N*(1520)0 --> p + pi-
  0.55/3.,       // id=65 N*(1520)0 --> n + pi0
  0.20/6.,       // id=65 N*(1520)0 --> Delta+ + pi-
  0.20/3.,       // id=65 N*(1520)0 --> Delta0 + pi0
  0.20/2.,       // id=65 N*(1520)0 --> Delta- + pi+
  0.20*2./3.,    // id=65 N*(1520)0 --> p + rho-
  0.20/3.,       // id=65 N*(1520)0 --> n + rho0
  0.0449*2./3.,  // id=65 N*(1520)0 --> n + pi+ + pi-
  0.0449/3.,     // id=65 N*(1520)0 --> n + pi0 + pi0
  0.0051,        // id=65 N*(1520)0 --> n + photon
  // NS110, 13 channels
  0.46/3.,       // id=66 N*(1535)0 --> p + pi-
  0.46*2./3.,    // id=66 N*(1535)0 --> n + pi0
  0.3875,        // id=66 N*(1535)0 --> n + eta
  0.01/2.,       // id=66 N*(1535)0 --> Delta+ + pi-
  0.01/3.,       // id=66 N*(1535)0 --> Delta0 + pi0
  0.01/6.,       // id=66 N*(1535)0 --> Delta- + pi+
  0.04/3.,       // id=66 N*(1535)0 --> p + rho-
  0.04*2/3.,     // id=66 N*(1535)0 --> n + rho0
  0.03*2./3.,    // id=66 N*(1535)0 --> n + pi+ + pi-
  0.03/3.,       // id=66 N*(1535)0 --> n + pi0 + pi0
  0.07/3.,       // id=66 N*(1535)0 --> N*(1440)+ + pi-
  0.07*2./3.,    // id=66 N*(1535)0 --> N*(1440)0 + pi0
  0.0025,        // id=66 N*(1535)0 --> n + photon
  // J/Psi, 8 channels
  0.0602,        // id=67 J/Psi --> e+ + e-
  0.0601,        // id=67 J/Psi --> mu+ + mu-
  0.0088,        // id=67 J/Psi --> dilepton + photon
  0.0337,        // id=67 J/Psi --> 2pi+ + 2pi- + pi0
  0.0290,        // id=67 J/Psi --> 3pi+ + 3pi- + pi0
  0.0150,        // id=67 J/Psi --> pi+ + pi- + pi0
  0.0120,        // id=67 J/Psi --> K+ + K- + pi+ + pi- + pi0
  0.7812,        // id=67 J/Psi --> junk
  // Psi' = Psi(2S), 8 channels
  0.0073,        // id=68 Psi'  --> e+e-
  0.0070,        // id=68 Psi'  --> mu+mu-
  0.31,          // id=68 Psi'  --> J/Psi pi+ pi-
  0.182,         // id=68 Psi'  --> J/Psi pi0 pi0
  0.027,         // id=68 Psi'  --> J/Psi eta
  0.0035,        // id=68 Psi'  --> 2pi+ + 2pi- + pi0
  0.0030,        // id=68 Psi'  --> 3pi+ + 3pi- + pi0
  0.4554,        // id=68 Psi'  --> junk
  // pn, 2 channels
  0.999,        // id=69 pn  --> dilepton + p + n
  0.001         // id=69 pn  --> photon + p + n
};

// Decay-mode nomenclature:
// id1 + id2*1000 + .. + idn*1000^(n-1), or id of non-virtual-photon if Dalitz
const char *PStdData::MODE[mnmodes]={
  // mu+, 1 channel
  "4004002",        // id=5  mu+ --> e+ + neutrino + neutrino
  // mu-, 1 channel
  "4004003",        // id=6  mu- --> e- + neutrino + neutrino
  // pi0, 2 channels
  "1001",           // id=7  pi0 --> photon + photon
  "1051",           // id=7  pi0 --> dilepton + photon (Dalitz)
  // pi+, 1 channel
  "4005",           // id=8  pi+ --> mu+ + neutrino
  // pi-, 1 channel
  "4006",           // id=9  pi- --> mu- + neutrino
  // K0L, 6 channels
  "7007007",        // id=10 K0 long --> pi0 + pi0 + pi0
  "7009008",        // id=10 K0 long --> pi+ + pi- + pi0
  "4006008",        // id=10 K0 long --> pi+ + mu- + neutrino
  "4005009",        // id=10 K0 long --> pi- + mu+ + neutrino
  "4003008",        // id=10 K0 long --> pi+ + e- + neutrino
  "4002009",        // id=10 K0 long --> pi- + e+ + neutrino
  // K+, 6 channels
  "4005",           // id=11 K+ --> mu+ + neutrino
  "8007",           // id=11 K+ --> pi+ + pi0
  "9008008",        // id=11 K+ --> pi+ + pi+ + pi-
  "7007008",        // id=11 K+ --> pi+ + pi0 + pi0
  "4005007",        // id=11 K+ --> pi0 + mu+ + neutrino
  "4002007",        // id=11 K+ --> pi0 + e+ + neutrino
  // K-, 6 channels
  "4006",           // id=12 K- --> mu- + neutrino
  "7009",           // id=12 K- --> pi- + pi0
  "8009009",        // id=12 K- --> pi- + pi- + pi+
  "7007009",        // id=12 K- --> pi- + pi0 + pi0
  "4006007",        // id=12 K- --> pi0 + mu- + neutrino
  "4003007",        // id=12 K- --> pi0 + e- + neutrino
  // n, 1 channel
  "4003014",        // id=13 n --> p + e- + neutrino
  // K0S, 3 channels
  "9008",           // id=16 K0 short --> pi+ + pi-
  "7007",           // id=16 K0 short --> pi0 + pi0
  "1009008",        // id=16 K0 short --> pi+ + pi- + photon
  // eta, 8 channels
  "1001",           // id=17 eta --> photon + photon
  "7007007",        // id=17 eta --> pi0 + pi0 + pi0
  "8009007",        // id=17 eta --> pi+ + pi- + pi0
  "8009001",        // id=17 eta --> pi+ + pi- + photon
  "1051",           // id=17 eta --> dilepton + photon (Dalitz)
  "1050",           // id=17 eta --> photon + dimuon
  "3002",           // id=17 eta --> e+ + e-
  "5006",           // id=17 eta --> mu+ + mu-
  // Lambda, 3 channels
  "9014",           // id=18 Lambda --> p + pi-
  "7013",           // id=18 Lambda --> n + pi0
  "1013",           // id=18 Lambda --> n + photon
  // anti_n, 1 channel
  "4002015",        // id=25 anti_n --> anti_p + e+ + neutrino
  // D0, 4 channels
  "9014",           // id=34 Delta0 --> p + pi-
  "7013",           // id=34 Delta0 --> n + pi0
  "1013",           // id=34 Delta0 --> n + photon
  "13051",          // id=34 Delta0 --> dilepton + n (Dalitz)
  // D++, 1 channel
  "8014",           // id=35 Delta++ --> p + pi+
  // D+, 4 channels
  "7014",           // id=36 Delta+ --> p + pi0
  "8013",           // id=36 Delta+ --> n + pi+
  "1014",           // id=36 Delta+ --> p + photon
  "14051",          // id=36 Delta+ --> dilepton + p (Dalitz)
  // D-, 1 channel
  "9013",           // id=36 Delta- --> n + pi-
  // NP11+, 10 channels
  "7014",           // id=38 N*(1440)+ --> p + pi0
  "8013",           // id=38 N*(1440)+ --> n + pi+
  "9035",           // id=38 N*(1440)+ --> Delta++ + pi-
  "7036",           // id=38 N*(1440)+ --> Delta+ + pi0
  "8034",           // id=38 N*(1440)+ --> Delta0 + pi+
  "41014",          // id=38 N*(1440)+ --> p + rho0
  "42013",          // id=38 N*(1440)+ --> n + rho+
  "7007014",        // id=38 N*(1440)+ --> p + pi0 + pi0
  "9008014",        // id=38 N*(1440)+ --> p + pi+ + pi-
  "1014",           // id=38 N*(1440)+ --> p + photon
  // ND13+, 10 channels
  "7014",           // id=39 N*(1520)+ --> p + pi0
  "8013",           // id=39 N*(1520)+ --> n + pi+
  "9035",           // id=39 N*(1520)+ --> Delta++ + pi-
  "7036",           // id=39 N*(1520)+ --> Delta+ + pi0
  "8034",           // id=39 N*(1520)+ --> Delta0 + pi+
  "41014",          // id=39 N*(1520)+ --> p + rho0
  "42013",          // id=39 N*(1520)+ --> n + rho+
  "7007014",        // id=39 N*(1520)+ --> p + pi0 + pi0
  "9008014",        // id=39 N*(1520)+ --> p + pi+ + pi-
  "1014",           // id=39 N*(1520)+ --> p + photon
  // NS11+, 13 channels
  "7014",           // id=40 N*(1535)+ --> p + pi0
  "8013",           // id=40 N*(1535)+ --> n + pi+
  "17014",          // id=40 N*(1535)+ --> p + eta
  "9035",           // id=40 N*(1535)+ --> Delta++ + pi-
  "7036",           // id=40 N*(1535)+ --> Delta+ + pi0
  "8034",           // id=40 N*(1535)+ --> Delta0 + pi+
  "41014",          // id=40 N*(1535)+ --> p + rho0
  "42013",          // id=40 N*(1535)+ --> n + rho+
  "7007014",        // id=40 N*(1535)+ --> p + pi0 + pi0
  "9008014",        // id=40 N*(1535)+ --> p + pi+ + pi-
  "7038",           // id=40 N*(1535)+ --> N*(1440)+ + pi0
  "8064",           // id=40 N*(1535)+ --> N*(1440)0 + pi+
  "1014",           // id=40 N*(1535)+ --> p + photon
  // rho0, 3 channels
  "9008",           // id=41 rho0 --> pi+ + pi-
  "3002",           // id=41 rho0 --> e+ + e-
  "5006",           // id=41 rho0 --> mu+ + mu-
  // rho+, 1 channel
  "7008",           // id=42 rho+ --> pi+ + pi0 
  // rho-, 1 channel
  "7009",           // id=43 rho- --> pi- + pi0 
  // dimuon, 1 channel
  "6005",           // id=50  dimuon --> mu+ + mu-
  // dilepton, 1 channel
  "3002",           // id=51 dilepton --> e+ + e-
  // w, 7 channels
  "7009008",        // id=52 omega --> pi+ + pi- + pi0
  "1007",           // id=52 omega --> pi0 + photon
  "9008",           // id=52 omega --> pi+ + pi-
  "7051",           // id=52 omega --> dilepton + pi0 (Dalitz)
  "7050",           // id=52 omega --> dimuon + pi0
  "3002",           // id=52 omega --> e+ + e-
  "5006",           // id=52 omega --> mu+ + mu-
  // eta', 7 channels
  "9008017",        // id=53 eta' --> eta + pi- + pi+ 
  "1041",           // id=53 eta' --> rho0 + photon
  "7007017",        // id=53 eta' --> eta + pi0 + pi0
  "1052",           // id=53 eta' --> omega + photon
  "1001",           // id=53 eta' --> photon + photon
  "7007007",        // id=53 eta' --> pi0 + pi0 + pi0
  "1050",           // id=53 eta' --> dimuon + photon
  // sigma, 3 channels
  "9008",           // id=54 sigma --> pi+ + pi-
  "3002",           // id=54 sigma --> e+ + e-
  "5006",           // id=54 sigma --> mu+ + mu-  
  // phi, 8 channels
  "12011",          // id=55 phi --> K+ + K-
  "16010",          // id=55 phi --> K0L + K0S
  "7009008",        // id=55 phi --> pi+ + pi- + pi0
  "1017",           // id=55 phi --> eta + photon
  "1007",           // id=55 phi --> pi0 + photon
  "3002",           // id=55 phi --> e+ + e-
  "5006",           // id=55 phi --> mu+ + mu-
  "17051",          // id=55 phi --> dilepton + eta
  // DP330, 10 channels
  "9014",           // id=56 Delta(1600)0 --> p + pi-
  "7013",           // id=56 Delta(1600)0 --> n + pi0
  "9036",           // id=56 Delta(1600)0 --> Delta+ + pi-
  "7034",           // id=56 Delta(1600)0 --> Delta0 + pi0
  "8037",           // id=56 Delta(1600)0 --> Delta- + pi+
  "43014",          // id=56 Delta(1600)0 --> p + rho-
  "41013",          // id=56 Delta(1600)0 --> n + rho0
  "43038",          // id=56 Delta(1600)0 --> N(1440)+ + rho-
  "41064",          // id=56 Delta(1600)0 --> N(1440)0 + rho0
  "1013",           // id=56 Delta(1600)0 --> n + photon
  // DP33++, 5 channels
  "8014",           // id=57 Delta(1600)++ --> p + pi+
  "7035",           // id=57 Delta(1600)++ --> Delta++ + pi0
  "8036",           // id=57 Delta(1600)++ --> Delta+ + pi+
  "42014",          // id=57 Delta(1600)++ --> p + rho+
  "42038",          // id=57 Delta(1600)++ --> N(1440)+ + rho+
  // DP33+, 10 channels
  "7014",           // id=58 Delta(1600)+ --> p + pi0
  "8013",           // id=58 Delta(1600)+ --> n + pi+
  "9035",           // id=58 Delta(1600)+ --> Delta++ + pi-
  "7036",           // id=58 Delta(1600)+ --> Delta+ + pi0
  "8034",           // id=58 Delta(1600)+ --> Delta0 + pi+
  "41014",          // id=58 Delta(1600)+ --> p + rho0
  "42013",          // id=58 Delta(1600)+ --> n + rho+
  "41038",          // id=58 Delta(1600)+ --> N(1440)+ + rho0
  "42064",          // id=58 Delta(1600)+ --> N(1440)0 + rho+
  "1014",           // id=58 Delta(1600)+ --> p + photon
  // DP33-, 5 channels
  "9013",           // id=59 Delta(1600)- --> n + pi-
  "9034",           // id=59 Delta(1600)- --> Delta0 + pi-
  "7037",           // id=59 Delta(1600)- --> Delta- + pi0
  "43013",          // id=59 Delta(1600)- --> n + rho-
  "43064",          // id=59 Delta(1600)- --> N(1440)0 + rho-
  // DS310, 8 channels
  "9014",           // id=60 Delta(1620)0 --> p + pi-
  "7013",           // id=60 Delta(1620)0 --> n + pi0
  "9036",           // id=60 Delta(1620)0 --> Delta+ + pi-
  "7034",           // id=60 Delta(1620)0 --> Delta0 + pi0
  "8037",           // id=60 Delta(1620)0 --> Delta- + pi+
  "43014",          // id=60 Delta(1620)0 --> p + rho-
  "41013",          // id=60 Delta(1620)0 --> n + rho0
  "1013",           // id=60 Delta(1620)0 --> n + photon
  // DS31++, 4 channels
  "8014",           // id=61 Delta(1620)++ --> p + pi+
  "7035",           // id=61 Delta(1620)++ --> Delta++ + pi0
  "8036",           // id=61 Delta(1620)++ --> Delta+ + pi+
  "42014",          // id=61 Delta(1620)++ --> p + rho+
  // DS31+, 8 channels
  "7014",           // id=62 Delta(1620)+ --> p + pi0
  "8013",           // id=62 Delta(1620)+ --> n + pi+
  "9035",           // id=62 Delta(1620)+ --> Delta++ + pi-
  "7036",           // id=62 Delta(1620)+ --> Delta+ + pi0
  "8034",           // id=62 Delta(1620)+ --> Delta0 + pi+
  "41014",          // id=62 Delta(1620)+ --> p + rho0
  "42013",          // id=62 Delta(1620)+ --> n + rho+
  "1014",           // id=62 Delta(1620)+ --> p + photon
  // DS31-, 4 channels
  "9013",           // id=63 Delta(1620)- --> n + pi-
  "9034",           // id=63 Delta(1620)- --> Delta0 + pi-
  "7037",           // id=63 Delta(1620)- --> Delta- + pi0
  "43013",          // id=63 Delta(1620)- --> n + rho-
  // NP110, 10 channels
  "9014",           // id=64 N*(1440)0 --> p + pi-
  "7013",           // id=64 N*(1440)0 --> n + pi0
  "9036",           // id=64 N*(1440)0 --> Delta+ + pi-
  "7034",           // id=64 N*(1440)0 --> Delta0 + pi0
  "8037",           // id=64 N*(1440)0 --> Delta- + pi+
  "43014",          // id=64 N*(1440)0 --> p + rho-
  "41013",          // id=64 N*(1440)0 --> n + rho0
  "9008013",        // id=64 N*(1440)0 --> n + pi+ + pi-
  "7007013",        // id=64 N*(1440)0 --> n + pi0 + pi0
  "1013",           // id=64 N*(1440)0 --> n + photon
  // ND130, 10 channels
  "9014",           // id=65 N*(1520)0 --> p + pi-
  "7013",           // id=65 N*(1520)0 --> n + pi0
  "9036",           // id=65 N*(1520)0 --> Delta+ + pi-
  "7034",           // id=65 N*(1520)0 --> Delta0 + pi0
  "8037",           // id=65 N*(1520)0 --> Delta- + pi+
  "43014",          // id=65 N*(1520)0 --> p + rho-
  "41013",          // id=65 N*(1520)0 --> n + rho0
  "9008013",        // id=65 N*(1520)0 --> n + pi+ + pi-
  "7007013",        // id=65 N*(1520)0 --> n + pi0 + pi0
  "1013",           // id=65 N*(1520)0 --> n + photon
  // NS110, 13 channels
  "9014",           // id=66 N*(1535)0 --> p + pi-
  "7013",           // id=66 N*(1535)0 --> n + pi0
  "17013",          // id=66 N*(1535)0 --> n + eta
  "9036",           // id=66 N*(1535)0 --> Delta+ + pi-
  "7034",           // id=66 N*(1535)0 --> Delta0 + pi0
  "8037",           // id=66 N*(1535)0 --> Delta- + pi+
  "43014",          // id=66 N*(1535)0 --> p + rho-
  "41013",          // id=66 N*(1535)0 --> n + rho0
  "9008013",        // id=66 N*(1535)0 --> n + pi+ + pi-
  "7007013",        // id=66 N*(1535)0 --> n + pi0 + pi0
  "9038",           // id=66 N*(1535)0 --> N*(1440)+ + pi-
  "7064",           // id=66 N*(1535)0 --> N*(1440)0 + pi0
  "1013",           // id=66 N*(1535)0 --> n + photon
  // J/Psi, 8 channels
  "2003",           // J/Psi --> e+ + e-
  "5006",           // J/Psi --> mu+ + mu-
  "1051",           // J/Psi --> dilepton + photon
  "8009008009007",  // J/Psi --> 2pi+ + 2pi- + pi0
  "8009008009008009007",  // J/Psi --> 3pi+ + 3pi- + pi0
  "8009007",        // J/Psi --> pi+ + pi- + pi0
  "12011008009007", // J/Psi --> K+ + K- + pi+ + pi- + pi0
  "4004",           // J/Psi --> 2 neutrinos  (this is a placeholder)
  // Psi' = Psi(2S), 8 channels
  "2003",           // id=68 Psi'  --> e+e-
  "5006",           // id=68 Psi'  --> mu+mu-
  "67008009",       // id=68 Psi'  --> J/Psi pi+ pi-
  "67007007",       // id=68 Psi'  --> J/Psi pi0 pi0
  "67017",          // id=68 Psi'  --> J/Psi eta
  "8009008009007",  // id=68 Psi'  --> 2pi+ + 2pi- + pi0
  "8009008009008009007",// id=68 Psi'  --> 3pi+ + 3pi- + pi0
  "4004",           // id=68 Psi'  --> junk (=2 neutrinos as a placeholder)
  // pn, 2 channels
  "51013014",       // id=69 pn  --> dilepton + p + n
  "1013014"         // id=69 pn  --> photon + p + n
};

// Decay-mode text description
const char *PStdData::DESCRIPTION[mnmodes]={
  "mu+ --> e+ + neutrino + neutrino",
  "mu- --> e- + neutrino + neutrino",
  "pi0 --> photon + photon",
  "pi0 --> dilepton + photon (Dalitz)",
  "pi+ --> mu+ + neutrino",
  "pi- --> mu- + neutrino",
  "K0 long --> pi0 + pi0 + pi0",
  "K0 long --> pi+ + pi- + pi0",
  "K0 long --> pi+ + mu- + neutrino",
  "K0 long --> pi- + mu+ + neutrino",
  "K0 long --> pi+ + e- + neutrino",
  "K0 long --> pi- + e+ + neutrino",
  "K+ --> mu+ + neutrino",
  "K+ --> pi+ + pi0",
  "K+ --> pi+ + pi+ + pi-",
  "K+ --> pi+ + pi0 + pi0",
  "K+ --> pi0 + mu+ + neutrino",
  "K+ --> pi0 + e+ + neutrino",
  "K- --> mu- + neutrino",
  "K- --> pi- + pi0",
  "K- --> pi- + pi- + pi+",
  "K- --> pi- + pi0 + pi0",
  "K- --> pi0 + mu- + neutrino",
  "K- --> pi0 + e- + neutrino",
  "n --> p + e- + neutrino",
  "K0 short --> pi+ + pi-",
  "K0 short --> pi0 + pi0",
  "K0 short --> pi+ + pi- + photon",
  "eta --> photon + photon",
  "eta --> pi0 + pi0 + pi0",
  "eta --> pi+ + pi- + pi0",
  "eta --> pi+ + pi- + photon",
  "eta --> dilepton + photon (Dalitz)",
  "eta --> dimuon + photon",
  "eta --> e+ + e-",
  "eta --> mu+ + mu-",
  "Lambda --> p + pi-",
  "Lambda --> n + pi0",
  "Lambda --> n + photon",
  "anti_n --> anti_p + e+ + neutrino",
  "Delta0 --> p + pi-",
  "Delta0 --> n + pi0",
  "Delta0 --> n + photon",
  "Delta0 --> dilepton + n (Dalitz)",
  "Delta++ --> p + pi+",
  "Delta+ --> p + pi0",
  "Delta+ --> n + pi+",
  "Delta+ --> p + photon",
  "Delta+ --> dilepton + p (Dalitz)",
  "Delta- --> n + pi-",
  "N*(1440)+ --> p + pi0",
  "N*(1440)+ --> n + pi+",
  "N*(1440)+ --> Delta++ + pi-",
  "N*(1440)+ --> Delta+ + pi0",
  "N*(1440)+ --> Delta0 + pi+",
  "N*(1440)+ --> p + rho0",
  "N*(1440)+ --> n + rho+",
  "N*(1440)+ --> p + pi0 + pi0",
  "N*(1440)+ --> p + pi+ + pi-",
  "N*(1440)+ --> p + photon",
  "N*(1520)+ --> p + pi0",
  "N*(1520)+ --> n + pi+",
  "N*(1520)+ --> Delta++ + pi-",
  "N*(1520)+ --> Delta+ + pi0",
  "N*(1520)+ --> Delta0 + pi+",
  "N*(1520)+ --> p + rho0",
  "N*(1520)+ --> n + rho+",
  "N*(1520)+ --> p + pi0 + pi0",
  "N*(1520)+ --> p + pi+ + pi-",
  "N*(1520)+ --> p + photon",
  "N*(1535)+ --> p + pi0",
  "N*(1535)+ --> n + pi+",
  "N*(1535)+ --> p + eta",
  "N*(1535)+ --> Delta++ + pi-",
  "N*(1535)+ --> Delta+ + pi0",
  "N*(1535)+ --> Delta0 + pi+",
  "N*(1535)+ --> p + rho0",
  "N*(1535)+ --> n + rho+",
  "N*(1535)+ --> p + pi0 + pi0",
  "N*(1535)+ --> p + pi+ + pi-",
  "N*(1535)+ --> N*(1440)+ + pi0",
  "N*(1535)+ --> N*(1440)0 + pi+",
  "N*(1535)+ --> p + photon",
  "rho0 --> pi+ + pi-",
  "rho0 --> e+ + e-",
  "rho0 --> mu+ + mu-",
  "rho+ --> pi+ + pi0",
  "rho- --> pi- + pi0",
  "dimuon --> mu+ + mu-",
  "dilepton --> e+ + e-",
  "omega --> pi+ + pi- + pi0",
  "omega --> pi0 + photon",
  "omega --> pi+ + pi-",
  "omega --> dilepton + pi0 (Dalitz)",
  "omega --> dimuon + pi0",
  "omega --> e+ + e-",
  "omega --> mu+ + mu-",
  "eta' --> eta + pi- + pi+ ",
  "eta' --> rho0 + photon",
  "eta' --> eta + pi0 + pi0",
  "eta' --> omega + photon",
  "eta' --> photon + photon",
  "eta' --> pi0 + pi0 + pi0",
  "eta' --> dimuon + photon",
  "sigma --> pi+ + pi-",
  "sigma --> e+ + e-",
  "sigma --> mu+ + mu-",
  "phi --> K+ + K-",
  "phi --> K0L + K0S",
  "phi --> pi+ + pi- + pi0",
  "phi --> eta + photon",
  "phi --> pi0 + photon",
  "phi --> e+ + e-",
  "phi --> mu+ + mu-",
  "phi --> eta + dilepton",
  "Delta(1600)0 --> p + pi-",
  "Delta(1600)0 --> n + pi0",
  "Delta(1600)0 --> Delta+ + pi-",
  "Delta(1600)0 --> Delta0 + pi0",
  "Delta(1600)0 --> Delta- + pi+",
  "Delta(1600)0 --> p + rho-",
  "Delta(1600)0 --> n + rho0",
  "Delta(1600)0 --> N(1440)+ + rho-",
  "Delta(1600)0 --> N(1440)0 + rho0",
  "Delta(1600)0 --> n + photon",
  "Delta(1600)++ --> p + pi+",
  "Delta(1600)++ --> Delta++ + pi0",
  "Delta(1600)++ --> Delta+ + pi+",
  "Delta(1600)++ --> p + rho+",
  "Delta(1600)++ --> N(1440)+ + rho+",
  "Delta(1600)+ --> p + pi0",
  "Delta(1600)+ --> n + pi+",
  "Delta(1600)+ --> Delta++ + pi-",
  "Delta(1600)+ --> Delta+ + pi0",
  "Delta(1600)+ --> Delta0 + pi+",
  "Delta(1600)+ --> p + rho0",
  "Delta(1600)+ --> n + rho+",
  "Delta(1600)+ --> N(1440)+ + rho0",
  "Delta(1600)+ --> N(1440)0 + rho+",
  "Delta(1600)+ --> p + photon",
  "Delta(1600)- --> n + pi-",
  "Delta(1600)- --> Delta0 + pi-",
  "Delta(1600)- --> Delta- + pi0",
  "Delta(1600)- --> n + rho-",
  "Delta(1600)- --> N(1440)0 + rho-",
  "Delta(1620)0 --> p + pi-",
  "Delta(1620)0 --> n + pi0",
  "Delta(1620)0 --> Delta+ + pi-",
  "Delta(1620)0 --> Delta0 + pi0",
  "Delta(1620)0 --> Delta- + pi+",
  "Delta(1620)0 --> p + rho-",
  "Delta(1620)0 --> n + rho0",
  "Delta(1620)0 --> n + photon",
  "Delta(1620)++ --> p + pi+",
  "Delta(1620)++ --> Delta++ + pi0",
  "Delta(1620)++ --> Delta+ + pi+",
  "Delta(1620)++ --> p + rho+",
  "Delta(1620)+ --> p + pi0",
  "Delta(1620)+ --> n + pi+",
  "Delta(1620)+ --> Delta++ + pi-",
  "Delta(1620)+ --> Delta+ + pi0",
  "Delta(1620)+ --> Delta0 + pi+",
  "Delta(1620)+ --> p + rho0",
  "Delta(1620)+ --> n + rho+",
  "Delta(1620)+ --> p + photon",
  "Delta(1620)- --> n + pi-",
  "Delta(1620)- --> Delta0 + pi-",
  "Delta(1620)- --> Delta- + pi0",
  "Delta(1620)- --> n + rho-",
  "N*(1440)0 --> p + pi-",
  "N*(1440)0 --> n + pi0",
  "N*(1440)0 --> Delta+ + pi-",
  "N*(1440)0 --> Delta0 + pi0",
  "N*(1440)0 --> Delta- + pi+",
  "N*(1440)0 --> p + rho-",
  "N*(1440)0 --> n + rho0",
  "N*(1440)0 --> n + pi+ + pi-",
  "N*(1440)0 --> n + pi0 + pi0",
  "N*(1440)0 --> n + photon",
  "N*(1520)0 --> p + pi-",
  "N*(1520)0 --> n + pi0",
  "N*(1520)0 --> Delta+ + pi-",
  "N*(1520)0 --> Delta0 + pi0",
  "N*(1520)0 --> Delta- + pi+",
  "N*(1520)0 --> p + rho-",
  "N*(1520)0 --> n + rho0",
  "N*(1520)0 --> n + pi+ + pi-",
  "N*(1520)0 --> n + pi0 + pi0",
  "N*(1520)0 --> n + photon",
  "N*(1535)0 --> p + pi-",
  "N*(1535)0 --> n + pi0",
  "N*(1535)0 --> n + eta",
  "N*(1535)0 --> Delta+ + pi-",
  "N*(1535)0 --> Delta0 + pi0",
  "N*(1535)0 --> Delta- + pi+",
  "N*(1535)0 --> p + rho-",
  "N*(1535)0 --> n + rho0",
  "N*(1535)0 --> n + pi+ + pi-",
  "N*(1535)0 --> n + pi0 + pi0",
  "N*(1535)0 --> N*(1440)+ + pi-",
  "N*(1535)0 --> N*(1440)0 + pi0",
  "N*(1535)0 --> n + photon",
  "J/Psi --> e+ + e-",
  "J/Psi --> mu+ + mu-",
  "J/Psi --> dilepton + photon",
  "J/Psi --> 2pi+ + 2pi- + pi0",
  "J/Psi --> 3pi+ + 3pi- + pi0",
  "J/Psi --> pi+ + pi- + pi0",
  "J/Psi --> K+ + K- + pi+ + pi- + pi0",
  "J/Psi --> junk",
  "Psi'  --> e+ + e-",
  "Psi'  --> mu+ + mu-",
  "Psi'  --> J/Psi + pi+ + pi-",
  "Psi'  --> J/Psi + pi0 + pi0",
  "Psi'  --> J/Psi + eta",
  "Psi'  --> 2pi+ + 2pi- + pi0",
  "Psi'  --> 3pi+ + 3pi- + pi0",
  "Psi'  --> junk",
  "pn --> dilepton + p + n",
  "pn --> photon + p + n"
};



// static members
int PStdData::maxnumpar=mnpar;                       // number of particles
int PStdData::maxnummodes=mnmodes;                   // number of decay modes

char ** PStdData::PName=(char**)NAME;                // particle names
double * PStdData::PMass=(double*)MASS;              // masses
double * PStdData::PWidth=(double*)WIDTH;            // widths
int * PStdData::PMeson=(int*)MESON;                  // meson (1) or not (0)
int * PStdData::PBaryon=(int*)BARYON;                // baryon number
int * PStdData::PLepton=(int*)LEPTON;                // lepton number
int * PStdData::PCharge=(int*)CHARGE;                // charge
int * PStdData::PJ=(int*)SPIN;                       // 2 x J
int * PStdData::PParity=(int*)PARITY;                // Parity
int * PStdData::PI=(int*)ISPIN;                      // 2 x I
int * PStdData::PNModes=(int*)NMODES;                // number of decay modes
int * PStdData::Pkf=(int*)PYTHIAKF;                  // Pythia6 kf code
double * PStdData::PBR=(double*)BRR;                 // branching ratio
char ** PStdData::PMode=(char**)MODE;                // decay mode coded by product ids
char ** PStdData::PMDescription=(char**)DESCRIPTION; // text description of decay mode


const long double PStdData::hbar=6.582122e-25; // units of (GeV s)


ClassImp(PStdData)



