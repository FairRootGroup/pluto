// Author: M.A. Kagarlis
// Written: 31.01.99
// Revised: 17/10/2006   R.H.
// PData Class Header
// (copied in part to PStdData) IF

#ifndef _PSTDDATA_H_
#define _PSTDDATA_H_

#include "TROOT.h"
#include "TArrayI.h"
#include "TArrayD.h"
#include <iostream>

#include "PStaticData.h"

using namespace std;



class PStdData;
PStdData* makeStdData();
PStdData& fStdData();

class PStdData : public TObject {
  
  

    
 public:

  //constructor
  PStdData();
  //destructor
  ~PStdData();

  Bool_t fillDataBase(void); //Copies the static entries into the PDataBase


 private:

 int disable;

 void resetPosition() {
     // Resets the array PPosition
     
     int i, nm, position=0;
     for (i=0;i<maxnumpar;++i) {
	 nm=PNModes[i];              // number of decay modes of particle pid=i
	 if (!nm) PPosition[i]=-1;   // signifies stable particle
	 else {
	     PPosition[i]=position;    // index of 1st decay mode for current particle
	     position+=nm;             // reposition for the next particle
	 }
     }
 }

 int *PPosition;


  static int maxnumpar, maxnummodes, *Pkf, *PMeson, *PBaryon, *PLepton,
    *PCharge, *PJ, *PParity, *PI, *PNModes, *intcache,
    cachesize, save, nfiles;
  static double *PMass, *PWidth, *PBR, *dblcache, scale;
  static char **PName, **PMDescription, **PMode;


  // local storage area
  static int *pmes_tmp, *pbar_tmp, *plep_tmp, *pchar_tmp, *pspin_tmp,
      *pparity_tmp, *pispin_tmp, *pnmod_tmp, *id_tmp, *pkf_tmp;
  static double *pmass_tmp, *pwidth_tmp, *pbr_tmp;
  static char **pnam_tmp, **pmdescr_tmp, **pmod_tmp;
  
  static const char *MESSAGE[];
  static const char *NAME[];
  static double MASS[];
  static double WIDTH[];
 static const int PYTHIAKF[];
 static const int MESON[];
 static const int BARYON[];
 static const int LEPTON[];
 static const int CHARGE[];
 static const int SPIN[];
 static const int PARITY[];
 static const int ISPIN[];
 static const int NMODES[];
 static double BRR[];
 static const char *MODE[];
 static const char *DESCRIPTION[];

 static const long double hbar;

  ClassDef(PStdData,0) //Pluto Particle Std Data Class
};

#endif // _PSTDDATA_H_







