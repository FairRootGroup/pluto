// Author: Romain Holzmann
// Written: 22.08.01
// PDiLepton Class Header

#ifndef _PDiLepton_H_
#define _PDiLepton_H_

#include "PParticle.h"
#include "PUtils.h"

class PDiLepton: public PParticle {

 protected:

  Int_t prodId;    // product id
  Float_t m1;      // lower edge of mass range
  Float_t m2;      // upper edge of mass range
  Float_t pt1;     // lower edge of pt range
  Float_t pt2;     // upper edge of pt range
  Float_t y1;      // lower edge of rapidity range
  Float_t y2;      // upper edge of rapidity range

 public:
  PDiLepton(Float_t m1, Float_t m2, Float_t pt1, Float_t pt2, Float_t y1,
            Float_t y2);

  int getParticleId() {return prodId;}
  Float_t getM1()  { return m1; }
  Float_t getM2()  { return m1; }
  Float_t getPt1() { return pt1; }
  Float_t getPt2() { return pt2; }
  Float_t getY1()  { return y1; }
  Float_t getY2()  { return y2; }
  Int_t isDilepton() { return 1; }
  void samplePartCM(double &px, double &py, double &pz, double &E); 
  virtual void Print(const Option_t* delme=NULL) const;
  virtual ~PDiLepton() { }

 protected:

ClassDef(PDiLepton,1) // Pluto DiLepton generator Class
};
#endif // _PDiLepton_H_





































