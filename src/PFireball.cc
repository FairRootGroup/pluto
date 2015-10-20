////////////////////////////////////////////////////////
//  Fireball Class implementation file
//
//  The Fireball class sets up multi-nucleon quasi particles
//  that may subsequently decay thermally via standard decay modes.
//  A fireball with the desired properties, including the decay
//  channels, must be available in PData before invoking this
//  constructor.
//
//                    Author:  Romain Holzmann
//                    Written: 28.04.00
//                    Revised: 21.06.00 MK
//                    Revised: 26.06.2001 A.Toia
//                    Revised: 05.09.2005 R.H.
//                    Revised: 28.03.2006 R.H.  (use TF2 for E,M sampling)
//                    Revised: 02.02.2007 R.H.  (add E p**n exp(-E/T) behavior)
//                    Revised: 16.02.2007 R.H.  (add angle-dependent temperature)
//                    Revised: 21.02.2007 R.H.  (add sampling from histogram)
//
// Ref 1: Gosset et al. PRC 16 (1977) 629
//        (number of projectile participant nucleons for
//         given impact parameter b (in fm))
////////////////////////////////////////////////////////


#include <TError.h>
#include "PUtils.h"
#include "PThermal.h"
#include "PFireball.h"
#include "PDistributionManager.h"

ClassImp(PFireball)

PFireball:: PFireball(char* particle, float AGeV, float t1, float t2, float f,
            float b, float a2, float a4, float w1, float w2, int sp) :
            PParticle(particle) {
//
// Thermal source with temperature(s) T1 (T2), frac*f(T1) + (1-frac)*f(T2),
// optional blast, optional polar anisotropies (A2,A4), optional flow (v1,v2).
//
  prodId = makeStaticData()->GetParticleID(particle);
  makeDistributionManager()->LinkDB(); //Fill data base with life

  //get primary model
  model= makeDynamicData()->GetParticleModel(prodId);

  SetID(500+prodId);
  quasistable = 1;
  //  quasistable = (PData::UMass(prodId)-PData::LMass(prodId) < 0.01); // (quasi)stable particle
  if (makeStaticData()->GetParticleTotalWidth(prodId) > THERMAL_WIDTH) 
      quasistable = 0;
  // This must be consistent with the definition in PThermal

  T1 = t1 > 0. ? t1 : 0.05;
  T2 = t2 >= 0. ? t2 : 0.;
  if (f <= 0.) frac = 1.;
  else frac = f;
  if (b<0.) {  // neg. blast flags power of E p**n term
    if (T2 > T1) T2 = T1;  // T1 = T(0 deg) and T2 = T(90 deg) with T1>=T2
    power1 = frac;
    power2 = -b;
    A2 = a2;
    A4 = a4;
    blast = 0.;
    frac = 1.;
  } else {
    if (f>1.) frac = 1.;
    power1 = power2 = 0.;
    A2 = a2;
    A4 = a4;
    blast = b < 1.0 ? b : 0.0;
  }
  v1 = w1;
  v2 = w2;
  Ap = 40.;
  At = 40.;
  prob = 0.1;
  flag = 0;
  meanN = 0.;
  spect = sp;
  fA=NULL;

  fE=NULL;
  pfE=NULL;
  afE=NULL;

  fE1=NULL;
  fE2=NULL;
  fE3=NULL;
  fMt=NULL;
  fHisto=NULL;
  trueThermal = kTRUE;
  sample_option = 0; //no rotation by default
  mesh_option= 0;   //no mesh by default
  updateFunctions();
  part=NULL;
  mt_fac=0.;
  sig=0.;
  didx_old = -1;
  setToMidrapidity(AGeV);   // move source to midrapidity for Ebeam is AGeV
}

void PFireball::SetEpsilon(Double_t e) {
    if (fE) fE->SetEpsilon(e);
    if (fE1) fE1->SetEpsilon(e);
    if (fE2) fE2->SetEpsilon(e);
    if (fE3) fE3->SetEpsilon(e);
}

void PFireball::SetNpx(Int_t my_npx) {
    npx = my_npx;
    if (fE) fE->SetNpx(npx);
    if (fE1) fE1->SetNpx(npx);
    if (fE2) fE2->SetNpx(npx);
    if (fE3) fE3->SetNpx(npx);
}

void PFireball::SetNpy(Int_t my_npy) {
    npy = my_npy;
    if (fE) fE->SetNpy(npy);
    if (fE1) fE1->SetNpy(npy);
    if (fE2) fE2->SetNpy(npy);
    if (fE3) fE3->SetNpy(npy);
}

void PFireball::sampleECM(Double_t &E, Double_t &M, int didx) {
    if (didx != didx_old) {
	//cout << didx << ":::" << didx_old << endl;
	didx_old = didx; //prevent random initialization each time
	if (didx<0) 
	    fE->SetParameter(5,-1.1);
	else 
	    fE->SetParameter(5,(double) didx);
    }
    if (!mesh_option) {
	fE->MakeIntegral();
	fE->GetRandom2(E,M);
    } else {
	if (!afE) {
	    //add an adaptive mesh
	    char *modname = new char[120];
	    char *desname = new char[120];
	    const char *prodname=makeStaticData()->GetParticleName(prodId);
	    sprintf(modname,"%s_fireball@%s/thermal",prodname,prodname);
	    sprintf(desname,"Thermal model for %s",prodname);
	    Info("sampleECM","Creating mesh for '%s'. This can take some time....",desname);
	    pfE = new PFunction(modname,desname,-1);
	    pfE->SetFunction(fE);
	    afE = new PAdaptiveMeshN (0x3, 2,pfE , 0.022);
	    afE->SetRange(0,fE->GetXmin(),fE->GetXmax());
	    afE->SetRange(1,fE->GetYmin(),fE->GetYmax());
	    afE->SetThreshold(1.5,0.0001);
	    afE->SetMCPoints(100);
	    afE->Divide(5,2);
	    //afE->Divide(5,0);
	    Info("sampleECM","...done");
	}
	afE->GetRandom();
	E=afE->GetArrayValue(0);
	M=afE->GetArrayValue(1);
    }
    if (sample_option) {
	double E_tmp = (E*0.5+2*M)*0.5;
	double M_tmp = (2*M-E*0.5)*0.5;
	E=E_tmp;
	M=M_tmp;
    }
    //cout << "exit sampleECM" << endl;
}

void PFireball::updateFunctions(int id) {

  if (id) prodId=id;
  SetID(500+prodId);
  Double_t M = makeStaticData()->GetParticleMass(prodId);
  Double_t Mmin = PData::LMass(prodId);
  Double_t Mmax = PData::UMass(prodId);
  Double_t Emin = Mmin;
  Double_t Emax = Mmax + 10.*(T1+T2)/(1.-blast);
  npx = 200;
  npy = 1000;

  if (quasistable) {
    Mmin = 0.999*makeStaticData()->GetParticleMass(prodId);  // mass axis is dummy
    Mmax = 1.001*makeStaticData()->GetParticleMass(prodId);
    Emin = M;
    Emax = Mmax + 20.*(T1+T2)/(1.-blast);
    npx=500;
    npy = 5;  // TF2 needs at least 5 points/axis
  }
 
  if (!fE) {
      //cout << Emin << ":" << Emax<< ":" <<Mmin<< ":" <<Mmax << endl;
    fE = new PF2("d2NdEdM",PThermal::d2NdEdM,Emin,Emax,Mmin,Mmax,7); // energy formula
    gROOT->GetListOfFunctions()->Remove(fE);
    fE->SetNpx(npx);
    fE->SetNpy(npy);

  } else fE->SetRange(Emin, Mmin, Emax, Mmax);
  fE->SetParameters(T1,T2,frac,blast,float(prodId),-1.1,sample_option); //option=0: normal operation 
  fE->Eval(0,0);
  //fE->Draw("surf");

  

  if (!fE1) {
    fE1 = new PF2("d2NdEdM1",PThermal::d2NdEdM1,Emin,Emax,Mmin,Mmax,6); // energy formula
    gROOT->GetListOfFunctions()->Remove(fE1);
    fE1->SetNpx(npx);
    fE1->SetNpy(npy);
  } else fE1->SetRange(Emin, Mmin, Emax, Mmax);
  fE1->SetParameters(T1,T2,frac,blast,float(prodId),-1.1);
  //fE1->Draw("same");

  if (!fE2) {
    fE2 = new PF2("d2NdEdM2",PThermal::d2NdEdM2,Emin,Emax,Mmin,Mmax,6); // energy formula
    gROOT->GetListOfFunctions()->Remove(fE2);
    fE2->SetNpx(npx);
    fE2->SetNpy(npy);
  } else fE2->SetRange(Emin, Mmin, Emax, Mmax);
  fE2->SetParameters(T1,T2,frac,blast,float(prodId),-1.1);
  //fE2->Draw("same");

  if (!fE3) {
    fE3 = new PF2("d2NdEdTheta",PThermal::d2NdEdTheta,Emin,Emax,0.,TMath::Pi(),7);
    gROOT->GetListOfFunctions()->Remove(fE3);
    fE3->SetNpx(npx);
    fE3->SetNpy(180);
  } else fE3->SetRange(Emin, 0., Emax,TMath::Pi());
  fE3->SetParameters(T1,T2,power1,power2,makeStaticData()->GetParticleMass(prodId),A2,A4);
  //fE3->Draw("same");

  if (!fA) {
    fA = new TF1("dNdTheta",PThermal::dNdTheta,0.,TMath::Pi(),2);  // polar formula
    gROOT->GetListOfFunctions()->Remove(fA);
    fA->SetNpx(180);
  }
  fA->SetParameters(A2,A4);

  if (!fMt) {
    fMt = new TF1("dNdMt",PThermal::dNdMt,M,Emax,4);               // mt formula
    gROOT->GetListOfFunctions()->Remove(fMt);
    fMt->SetNpx(500);
  } else fMt->SetRange(M, Emax);
  fMt->SetParameters(M,T1,T2,frac);

}

void PFireball::setToMidrapidity(float agev) {  // boost thermal source to
  // the correct rapidity with agev the energy/nucleon of the projectile (in GeV/u):
  // 0 for target-like spectators, beam rapidity for projectile-like spectators,
  // midrapidity for participants
  if (agev <= 0.0) return;
  if (spect==0) {   
    double bx = 0.;
    double by = 0.;
    double bz = sqrt(agev/(agev+2.*0.9315));
    Boost(bx,by,bz);
  }
 else if (spect==1) return;
 else if (spect==2) {
    double bx = 0.; 
    double by = 0.;
    double bz = sqrt(agev*(agev+2*0.9315))/(agev+0.9315); 
    Boost(bx,by,bz);
  }
}

void PFireball:: Print(const Option_t* delme) const {
  if (fHisto!=NULL) {
    printf("%s d2N/dpdTheta sampled from histogram: %s\n",delme,fHisto->GetTitle());
  }
  else if(power1==0) {
    printf("%s  Thermal %s with:\n%s  T1=%5.3f  frac=%5.3f  T2=%5.3f  blast=%5.3f \n",delme,
	 makeStaticData()->GetParticleName(prodId),delme, T1, frac, T2, blast);
  } else {
    printf("%s  Thermal %s with:\n%s  T1=%5.3f  pow1=%3.1f  T2=%5.3f  pow2=%3.1f \n",delme,
	 makeStaticData()->GetParticleName(prodId),delme, T1, power1, T2, power2);
  }
  printf("%s  A2=%5.2f  A4=%5.2f  v1=%5.2f  v2=%5.2f\n", delme,A2, A4, v1, v2);
  if(flag) printf("%s   Aproj=%4.0f  Atarg=%4.0f  Pprod=%8.2e\n", delme,Ap, At, prob); 
  if (spect==0) {
    printf("%s  participant\n",delme);
    printf("%s  rapidity=%f;  beta=%f;  weight=%f;  mult=%f\n", delme,Rapidity(), Beta(), W(), mult);
  }

  if (spect==1) {
    printf("%s  target-like spectator\n",delme);
    printf("%s  no Boost; rapidity=%f;  weight=%f;  mult=%f\n", delme,Rapidity(), W(), mult);
  }

  if (spect==2) {
    printf("%s  projectile-like spectator\n",delme);
    printf("%s  rapidity=%f;  beta=%f;  weight=%f;  mult=%f\n", delme,Rapidity(), Beta(), W(), mult);
  }

  if (flag) {
      printAverages();   // b sampling is on
  } else {
      printf("%s  meanN=%f\n",delme,meanN);
  }
}

void PFireball:: printAverages() const {
  Float_t b;
  Float_t avb=0.;
  Float_t Apart;
  Float_t avApart=0.;
  Int_t nSample = 10000;
  for (Int_t i=0;i<nSample;i++) {
    b = sampleB();
    avb += b;
    Apart = NparSmeared(Ap,At,b) + NparSmeared(At,Ap,b);    // av. nb. of participants
    avApart += Apart;
  }
  avb /= (Float_t)nSample;
  avApart /= (Float_t)nSample;
  Float_t avM=avApart*prob;
  
  printf(" impact parameter range: %5.2f - %5.2f fm\n",bmin,bmax);
  printf(" average b=%5.2f  average Apart=%5.1f  average M=%8.2e\n\n", avb, avApart, avM);
}

Double_t PFireball:: AvApart(Double_t ap, Double_t am, Double_t bm, Double_t bx) {
  Int_t nSample = 10000;
  Float_t b;
  Float_t Apart;
  Float_t avApart=0.;
  for (Int_t i=0;i<nSample;i++) {
    b = sampleB();
    Apart = Npar(Ap,At,b) + Npar(At,Ap,b);    // av. nb. of participants
    avApart += Apart;
  }
  return avApart/(Float_t)nSample;
}

void PFireball:: samplePartCM(double& px, double& py, double& pz, double& E, int didx) {
  // sample particle 4-momentum (px,py,pz,E)
   Double_t E1, E2;
   Double_t M, M1, M2;
   Double_t p;
   Double_t phi, phiPlane=*(makeStaticData()->GetBatchValue("_event_plane"));
   Int_t count = 0;

   Double_t theta;

   if (sig>0.) {     // do mt and y sampling separately
     Double_t mt = sampleMt();                 // dN/dMt = mt^2 * K1(mt/T)
     Double_t y = PUtils::sampleGaus(0.,sig);  // dN/dy = exp(-0.5*(y/sig)^2)
     if (model) { 
	 do {
	     model->SampleMass(&M);
	     count++;
	 } while(M>mt && count<100);

	 if (count>=100) Warning("samplePartCM","Model mass sampling failed (count>=100)");

     } else { //no model -> stable mass
	 M=makeStaticData()->GetParticleMass(prodId);
     }
     Double_t pt = sqrt(mt*mt-M*M);
     pz = mt*sinh(y);     
     E = mt*cosh(y);
     theta = atan(pt/pz);
     //phi = 2.*TMath::Pi()*PUtils::sampleFlat();  // --> Flat sampling, no v1,v2
     Double_t cost = cos(theta);
     Double_t sint = sin(theta);
     Double_t val, valmax;
     Int_t i = 0;
     do {           // sample phi now
	 i++;
	 phi = 2.*TMath::Pi()*PUtils::sampleFlat(); 
	 val = 1. + 2.*(v1*cos(phi)*cost + v2*cos(2.*phi)*sint);
	 valmax = 1. + 2.*(TMath::Abs(v1)*cost + TMath::Abs(v2)*sint);
     } while (valmax*PUtils::sampleFlat()>val  &&  i<100);
     
     px = pt*cos(phi+phiPlane);
     py = pt*sin(phi+phiPlane);

     return;
   }



   Int_t saveLevel = gErrorIgnoreLevel;
   //if (gErrorIgnoreLevel<=1000) gErrorIgnoreLevel = 1001; // suppress TF2 warnings

   if(fHisto != NULL) {  // if histo defined use it
     M = makeStaticData()->GetParticleMass(prodId);
     fHisto->GetRandom2(p,theta);  // sample thermal distribution from 2-d histogram
     E = sqrt(p*p+M*M);
     theta = 0.0174532925*theta;
   }
   else if (trueThermal) {
     theta = sampleThetaCM();
     sampleECM(E,M, didx); // sample true thermal distribution
     if (quasistable) M = makeStaticData()->GetParticleMass(prodId);   // reset to exact pole mass
     p = sqrt(TMath::Abs(E*E-M*M));
   }
   else if (power1==0) {
     theta = sampleThetaCM();
     sampleECM1(E1,M1, didx);  // sample thermal distribution (E*sqrt(p))
     sampleECM2(E2,M2, didx);  // sample thermal distribution (with E*p**3)
     Double_t x = 2.*theta/TMath::Pi();
     if (x>1.) x = 2. - x;
     if (x<0.) x = 0.;
     if (x>1.) x = 1.;
     //  x  = sqrt(x);
     if (x>PUtils::sampleFlat()) {E = E1; M = M1;}
     else {E = E2; M = M2;}   // cook up an angle dependency of E
     if (quasistable) M = makeStaticData()->GetParticleMass(prodId);   // reset to exact pole mass
     p = sqrt(TMath::Abs(E*E-M*M));
   } else {
       M = makeStaticData()->GetParticleMass(prodId);
     sampleECM3(E,theta);  // sample thermal distribution E p**n with T=T(th) and n=n(th)
     p = sqrt(TMath::Abs(E*E-M*M));
   }

   gErrorIgnoreLevel = saveLevel;

   Double_t cost = cos(theta);
   Double_t sint = sin(theta);

   Double_t val, valmax;
   Int_t i = 0;
   do {           // sample phi now
      i++;
      phi = 2.*TMath::Pi()*PUtils::sampleFlat(); 
      val = 1. + 2.*(v1*cos(phi)*cost + v2*cos(2.*phi)*sint);
      valmax = 1. + 2.*(TMath::Abs(v1)*cost + TMath::Abs(v2)*sint);
   } while (valmax*PUtils::sampleFlat()>val  &&  i<100);

   px = p*sint*cos(phi+phiPlane);
   py = p*sint*sin(phi+phiPlane);
   pz = p*cost;
   //cout << "exit samplePartCM" << endl;
}

float PFireball:: sampleB()  const {
   if(!flag) return 0.;
   if (bmin>0.1) return sqrt((bmax*bmax-bmin*bmin)*PUtils::sampleFlat() + bmin*bmin); // random b
   float b, test;
   do {
      b = (bmax+2.)*sqrt(PUtils::sampleFlat());
      test = 1./(1.+TMath::Exp((b-bmax)/0.5));        // nuclear density distribution
   } while(test<PUtils::sampleFlat());
   return b;
}

int PFireball:: sampleNProd() {
   return PUtils::samplePoisson(meanN);
}

int PFireball:: sampleNProd(float b) {
   if(!flag || b > 1.14*(pow((double)Ap,(double)0.3333)+pow((double)At,(double)0.3333))) return 0;
   float Apart = 0.;
   if (spect==0) Apart = NparSmeared(Ap,At,b) + NparSmeared(At,Ap,b); // av. nb. of participants
   else if (spect==1)  Apart = At - Npar(At,Ap,b); // av. nb. of target-like spectators
   else if (spect==2)  Apart = Ap - Npar(Ap,At,b); // av. nb. of projectile-like spectators

   float Nprod = prob*Apart;                       // av. nb. of produced part.
  
   if(prob < 0.25) nProd = PUtils::samplePoisson(Nprod);
                              // Poisson distributed nb. of product particles
   else nProd = PUtils::sampleBinomial((Int_t)(Apart+0.5),prob);
                              // use binomial to make sure that nProd <= Apart
   return nProd;
}

float PFireball:: Npar(float ap, float at, float b) {
//    number of projectile participant nucleons
//    for given impact parameter b (in fm) (Ref 1)

   float Rp = 1.14*pow((double)ap,(double)0.3333); // in fm
   float Rt = 1.14*pow((double)at,(double)0.3333);
   if(b > Rp+Rt) return 0;
   float beta = b/(Rp+Rt);
   float nu = Rp/(Rp+Rt);
   float mu = Rt/Rp;

   float F = 0.;
   if (nu > 0.5) {                  // Ap > At
     if (nu > 0.5*(1.0+beta)) {
       F = (1.0-pow(1.-mu*mu,1.5)) * sqrt(1.-(beta/nu)*(beta/nu));
     } else {
       F = 0.75*sqrt(1.-nu) * (1.-beta)/nu * (1.-beta)/nu
            - 0.125*( 3.*sqrt(1.-nu)/mu
            - (1.-pow(1.-mu*mu,1.5)) * sqrt(1.-(1.-mu)*(1.-mu))/(mu*mu*mu) )
                *pow((1.-beta)/nu,3.);
     }
   } else {                        // Ap < At
     if(nu > 0.5*(1.-beta)) {
       F = 0.75*sqrt(1.-nu) * (1.-beta)/nu * (1.-beta)/nu
               -0.125*(3.*sqrt(1.-nu)-1.) * pow((1.-beta)/nu,3.);
     } else {
       F = 1.;
     }
   }
   return F*ap;
}

float PFireball:: NparSmeared(float ap, float at, float b)  const {
//    number of projectile participant nucleons
//    for given impact parameter b (in fm) with smeared nuclear surface

   float Rp = 1.14*pow((double)ap,(double)0.3333)+1.; // + 1 fm tail
   float Rt = 1.14*pow((double)at,(double)0.3333)+1.; // + 1 fm tail
   if(b > Rp+Rt) return 0;
   float beta = b/(Rp+Rt);
   float nu = Rp/(Rp+Rt);
   float mu = Rt/Rp;

   float F = 0.;
   if (nu > 0.5) {                  // Ap > At
     if (nu > 0.5*(1.0+beta)) {
       F = (1.0-pow(1.-mu*mu,1.5)) * sqrt(1.-(beta/nu)*(beta/nu));
     } else {
       F = 0.75*sqrt(1.-nu) * (1.-beta)/nu * (1.-beta)/nu
            - 0.125*( 3.*sqrt(1.-nu)/mu
            - (1.-pow(1.-mu*mu,1.5)) * sqrt(1.-(1.-mu)*(1.-mu))/(mu*mu*mu) )
                *pow((1.-beta)/nu,3.);
     }
   } else {                        // Ap < At
     if(nu > 0.5*(1.-beta)) {
       F = 0.75*sqrt(1.-nu) * (1.-beta)/nu * (1.-beta)/nu
               -0.125*(3.*sqrt(1.-nu)-1.) * pow((1.-beta)/nu,3.);
     } else {
       F = 1.;
     }
   }
   return F*ap;
}

PChannel* PFireball::makeChannel(Int_t nMax, Float_t nMean) {
  //
  // set up a reaction channel for this thermal source
  //
   if (nMean>0.) setMeanN(nMean);

   part = new PParticle*[nMax+1];
   part[0] = this;
   for (Int_t i=1;i<=nMax;i++) {
     part[i] = new PParticle((char*)makeStaticData()->GetParticleName(prodId)); 
   }
   PChannel* chan = new PChannel(part,nMax,1);
   return chan;
}

Float_t PFireball::mtIntegral(Double_t m, Float_t T) { // integral of thermal
                                                       // mt distribution
   if (m<=0. || T<=0.) return 0.;
   return T*exp(-m/T)*sqrt(m)*(1.5*T+m)
            + 1.329340*pow((double)T,(double)2.5)*(1.-TMath::Erf(sqrt(m/T)));
}













