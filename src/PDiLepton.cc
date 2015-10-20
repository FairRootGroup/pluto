////////////////////////////////////////////////////////
//  Dilepton generator Class implementation file
//
//  The Dilepton class sets up an e+e- quasiparticle
//  that subsequently decays into e+ and e-.
//  Mass, pt and rapidity are sampled from a box distribution.
// 
//                    Author:  Romain Holzmann
//                    Written: 22.08.01
//
////////////////////////////////////////////////////////

#include "PDiLepton.h"

ClassImp(PDiLepton)

PDiLepton:: PDiLepton(Float_t M1, Float_t M2, Float_t Pt1, Float_t Pt2,
                      Float_t Y1, Float_t Y2) : PParticle("dilepton") {
//
// Dilepton source with m1<mass<m2, pt1<pt<pt2 and y1<rapidity<y2
//
  prodId = makeStaticData()->GetParticleID("dilepton");
  SetID(500+prodId);
  m1 = M1;
  m2 = M2;
  if (m1 > m2) { m2 = M1; m1 = M2; }
  if (m1 < 0.) m1= 0.;
  pt1 = Pt1;
  pt2 = Pt2;
  if (pt1 > pt2) { pt2 = Pt1; pt1 = Pt2; }
  if (pt1 < 0.) pt1 = 0.;
  y1 = Y1;
  y2 = Y2;
  if (y1 > y2) { y2 = Y1; y1 = Y2; }
}

void PDiLepton:: Print(const Option_t* delme) const{
  printf(" Dilepton with:\n %5.2f < M < %5.2f,  %5.2f < Pt < %5.2f and %5.2f < y < %5.2f\n",
	 m1,m2,pt1,pt2,y1,y2);
}

void PDiLepton:: samplePartCM(double& px, double& py, double& pz, double& E) {
  // sample particle 4-momentum (px,py,pz,E)
  Double_t m;
  Double_t y;
  Double_t pt;
  Double_t phi;

  do {  
    m = m1 + (m2-m1)*PUtils::sampleFlat();     // sample m, pt, y and phi uniformly
  } while (m < makeStaticData()->GetParticleMass(prodId));
  pt = pt1 + (pt2-pt1)*PUtils::sampleFlat();
  y = y1 + (y2-y1)*PUtils::sampleFlat();
  phi = 2.*TMath::Pi()*PUtils::sampleFlat();

  Double_t mt = sqrt(pt*pt+m*m);
  E = mt*cosh(y);
  px = pt*cos(phi);
  py = pt*sin(phi);
  pz = mt*sinh(y);
  return;
}















