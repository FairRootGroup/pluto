// Author: K.Tyminska
// Written: 20/11/2001
// Revised: 08/07/2004  R.H.
// 

#ifndef _PCROSS_H_
#define _PCROSS_H_
#include "TObject.h"
#include "TF1.h"
using namespace std;
#include <iostream>
#include <iomanip>
#include <cstdlib>

class PCross : public TObject {

 public:
  PCross();
  ~PCross(){;}
  static void print(int, Double_t bmin=0., Double_t bmax=100.);
  static void print(char *,Double_t bmin=0., Double_t bmax=100.);
  static void setSystem(int,int,int,int,Double_t Energy = 1.,Bool_t flag=0); 
  static void plot(int,Double_t,Double_t, Double_t bmin=0., 
                           Double_t bmax=100., const char * Opt = "L", Int_t col=1);
  static void plot(char *,Double_t,Double_t,Double_t bmin=0.,
                           Double_t bmax=100., const char * Opt = "L", Int_t col=1);
  static Double_t cross(int, Double_t bmin=0., Double_t bmax=100.);
  static Double_t cross(char *,Double_t bmin=0., Double_t bmax=100.);
  static Double_t calcT(Double_t);

 private:
  static Double_t Ebeam;   // in GeV/u
  static Double_t sqrts;   // in GeV
  static Bool_t sys;
  static Int_t  AP;
  static Int_t  ZP;
  static Int_t  AT;
  static Int_t  ZT;
  static Bool_t doMult;

  static Double_t Npart(int,int, Double_t);
  static Double_t ratiosignew(Double_t,Double_t); // T in MeV, mm = MeV 
  static Double_t calc(Double_t* , Double_t*);
  static Bool_t check(Int_t part);
 
ClassDef(PCross,0) // meson production cross sections in HIC
};
#endif // _PCROSS_H_



