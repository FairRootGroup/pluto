///////////////////////////////////////////////////////////////////////////////
//  Pluto cross section class
//  This utility class computes number of participants according to Gosset
//  and meson production probabilities for pi0, eta, omega and phi. 
//  pi0 and eta are averages of C+C and Ca+Ca TAPS (and KaoS) data,
//  omega and phi are mt-scaled from the eta averages.
//
//                             Author:  K.Tyminska
//                             Written: 20/12/2001 
//                             Revised: 08/07/2004   R.H.
//   Better calc. of source T  Revised: 04/09/2006   R.H.
//
// Ref 1: TAPS data is used (PRC 56 (1997) R2920)
// Ref 2: number of projectile participant nucleons for given impact
//        parameter b according to Gosset et al. PRC 16 (1977) 629
// Ref 3: chemical freeze-out temperature Cleymans et al. PRC 73 (2006) 034905
///////////////////////////////////////////////////////////////////////////////

#include "PCross.h"    
#include "PData.h"

Int_t PCross::ZP=0;
Int_t PCross::AP=0;
Int_t PCross::ZT=0;
Int_t PCross::AT=0;
Double_t PCross::Ebeam=0.0;
Double_t PCross::sqrts=0.0;
Bool_t PCross::sys=0;
Bool_t PCross::doMult=0;

ClassImp(PCross)

PCross::PCross() {sys=0;}

void PCross::print(int Part=0, Double_t blow, Double_t bup)
{  
  Double_t pmes, sigmes;

  if (check(Part)) {

//  empirical fireball temperature
     
    Double_t bmin = 0.;
    Double_t bmax = 1.14*(pow(AP,0.3333)+pow(AT,0.3333));   // in fm
    Double_t sigmax = 3.14159*bmax*bmax*10.;   // in mb
    Double_t T = 1000.*calcT(Ebeam);  // T in MeV

    if (bup>bmax) bup = bmax;
    bmin = blow;
  
    Double_t sigreac = 3.1416*(bup*bup - blow*blow)*10.;  // in mb
    Double_t frac = sigreac/sigmax*100;   // in %  

//  average Apart
    Double_t Apart = (AP*pow(AT,0.667) + AT*pow(AP,0.667))
                   /pow((pow(AP,0.333)+pow(AT,0.333)),2);

    int Np = AP - ZP;
    int Nt = AT - ZT;
    int N = Np + Nt;
    int Z = ZT + ZP;

//  ratio of pi-/pi0 in isobar model
    Double_t pim = (5.*N*(N-1.) + N*Z)/(Z*(Z-1.) + N*(N-1) + 4.*N*Z);
//  ratio of pi+/pi0 in isobar model
    Double_t pip = (5.*Z*(Z-1.) + N*Z)/(Z*(Z-1.) + N*(N-1) + 4.*N*Z);  

    Double_t bstep = (bup-blow)/100.;  // in fm
    Double_t sum1 = 0.;
    Double_t sum2 = 0.;

    for (Int_t i=1;i<101;i++){    
      Double_t b=blow + i*bstep;
      sum1 = sum1+b;
      sum2 = sum2 + b*(Npart(AP,AT,b)+Npart(AT,AP,b));  // Gosset function
    }
    Double_t NGosset = sum2/sum1;

    sigmes = cross(Part,blow,bup);   // compute cross section in mb
    pmes = sigmes/(sigreac*NGosset);

    cout << endl;
    cout <<"For reaction ("<<AP<<","<<ZP<<") + ("<<AT<<","<<ZT<<") at "<<Ebeam
         <<" AGeV"<<" (sqrt(s) ="<<sqrts<<") with b = "<<bmin<<"-"<<bmax<<" fm:"<<endl<<endl; 
    cout <<" <A>geom = "<<Apart<<" <N>Gosset = "<< NGosset<< endl;
    cout <<" sigR = "<< sigreac<<" mb  ("<<frac<<"% of total)"<<   endl;
    cout <<" Temp at "<< Ebeam<<" AGeV = "<<T<<" MeV"<< endl<<endl;

    if ((Part==7)||(Part==8)||(Part==9)){

      cout <<" pi-/pi0/pi+ = "<<pim<<":1:"<<pip<< endl;
      if (Part==7){
        cout <<" pi0    : Prob/part = "<<pmes<<"    Mult = "<<pmes*NGosset
                            <<"    Sig = "<<sigmes <<" mb" << endl;
      } else if (Part==8) {
        cout <<" pi+    : Prob/part = "<<pmes<<"    Mult = "<<pmes*NGosset
                         <<"    Sig = "<<sigmes <<" mb" << endl;
      } else {
        cout <<" pi-    : Prob/part = "<<pmes<<"    Mult = "<<pmes*NGosset
                         <<"    Sig = "<<sigmes <<" mb" << endl;
      }
    } else if (Part==17) {
      cout <<" eta    : Prob/part = "<<pmes<<"   Mult = "<<pmes*NGosset
	                 <<"    Sig = "<<sigmes <<" mb" << endl;
    } else {
      cout <<makeStaticData()->GetParticleName(Part)<<"    : Prob = "<<pmes<<" Mult = "<<pmes*NGosset
                         <<"  Sig = "<<sigmes<<" mb" << endl;
    }
  }
}


Double_t PCross::Npart(int AP,int AT,Double_t b)

//    returns number of projectile participant nucleons for given impact
//    parameter b according to Gosset et al. PRC 16 (1977) 629.
 {
  Double_t pom=0.0;
  Double_t Rp = 1.14*pow(AP,0.3333);
  Double_t Rt = 1.14*pow(AT,0.3333);
  Double_t beta=b/(Rp+Rt);
  Double_t nu = Rp/(Rp+Rt);
  Double_t mu = Rt/Rp;
  Double_t F;
  if (nu>0.5) {    //  Ap>At
    if (nu>0.5*(1.0+beta)) {
      F= (1.0-pow((1.0-pow(mu,2)),1.5))*sqrt(1.0-pow((beta/nu),2));      
    } else {
      F = 0.75*sqrt(1.0-nu)*pow(((1.0-beta)/nu),2) 
          - 0.125*(3.0*sqrt(1.0-nu)/mu - 
           ((1.0-pow((1.0-pow(mu,2)),1.5))*
           sqrt(1.0-pow((1.0-mu),2)))/pow(mu,3))*pow(((1.0-beta)/nu),3);
    }
  } else {        //  Ap<At 
    if (nu>0.5*(1.0-beta)) {
      F = 0.75*sqrt(1.0-nu)*pow(((1.0-beta)/nu),2) -
          0.125*(3.0*sqrt(1.0-nu) - 1.0)*pow(((1.0-beta)/nu),3);
    } else {
      F = 1.0;
    }
  }
  pom = F*AP;
  return pom;
}


Double_t PCross::ratiosignew(Double_t T,Double_t mm) // T in MeV, mm = MeV
  //    mt scaling applied to the eta and a meson of mass m(1)
  //    using integral[m,infinity] of dsigma/dmt = 1/T*mt**3/2*exp(-mt/T) 

{

  Double_t a=T*exp(-mm/T)*sqrt(mm)*(1.5*T+mm) 
                  + 0.75*1.772454*pow(T,2.5)*(1.-TMath::Erf(sqrt(mm/T)));
  Double_t b= T*exp(-547./T)*sqrt(547.)*(1.5*T+547.)
                  + 0.75*1.772454*pow(T,2.5)*(1.-TMath::Erf(sqrt(547./T)));
  Double_t pom = a/b;
  return pom;
}


void PCross::setSystem(int lAp,int lZp,int lAt,int lZt,Double_t lEbeam,Bool_t flag) {
  AP = lAp;
  ZP = lZp;
  AT = lAt;
  ZT = lZt;
  Ebeam = lEbeam;
  sqrts = sqrt(2.*0.939*(Ebeam+2.*0.939));  // NN sqrt(s)
  doMult = flag;
  sys = 1;
}


void PCross::print(char * Part, Double_t blow, Double_t bup) {
  return PCross::print(makeStaticData()->GetParticleID(Part),blow,bup);
}


void PCross::plot(char * Part,Double_t Rangl,Double_t Rangu, 
		  Double_t blow, Double_t bup,const char * Opt, Int_t color) {

  PCross::plot(makeStaticData()->GetParticleID(Part), Rangl,Rangu,blow,bup,Opt,color);
}


void PCross::plot(int Part,Double_t Rangl,Double_t Rangu, 
		   Double_t blow, Double_t bup,const char *Opt, Int_t color) {

  if (check(Part)) {
    TF1 *fn1 = new TF1("fn1",PCross::calc,Rangl,Rangu,7);
    fn1 -> SetParameters(Part,blow,bup,AT,ZT,AP,ZP);
    fn1 -> SetLineWidth(3);
    fn1 -> SetLineColor(color);    
    fn1 -> DrawCopy(Opt);
    delete fn1;
  }
  return;
}


Double_t PCross::calc(Double_t *x, Double_t *par) {
// return particle inclusive production cross section in mb
//
   Double_t Ene = x[0];
   int Part=(int)par[0];
   Double_t blow=par[1];
   Double_t bup=par[2];
   Int_t At=(Int_t)par[3]; // needed if more than 1 TF1 object is to be plotted
   Int_t Zt=(Int_t)par[4]; // cannot use static AP,ZP etc. because TF1s are
   Int_t Ap=(Int_t)par[5]; // recalculated on TCanvas::Update() ...
   Int_t Zp=(Int_t)par[6];

   Double_t T = 1000.*calcT(Ene); //  empirical meson temperature in MeV
   Double_t bmin = 0.;
   Double_t bmax = 1.14*(pow(Ap,0.3333)+pow(At,0.3333));  // in fm
   Double_t sigmax = 3.14159*bmax*bmax*10.;  // in mb 

   if (bup>bmax) bup = bmax;
   bmin = blow;
  
   Double_t sigreac = 3.1416*(bup*bup - blow*blow)*10.;  // in mb 
   if (doMult) sigreac = 1.;

   Double_t Apart = (Ap*pow(At,0.667) + At*pow(Ap,0.667))
                   /pow((pow(Ap,0.333)+pow(At,0.333)),2); // average Apart

   int Np = Ap - Zp;
   int Nt = At - Zt;
   int N = Np + Nt;
   int Z = Zt + Zp;

// pion ratios in isobar model
   Double_t pim = (5.*N*(N-1.) + N*Z)/(Z*(Z-1.) + N*(N-1) + 4.*N*Z);
   Double_t pip = (5.*Z*(Z-1.) + N*Z)/(Z*(Z-1.) + N*(N-1) + 4.*N*Z);  

   Double_t bstep = (bup-blow)/100.;  // in fm
   Double_t sum1 = 0.;
   Double_t sum2 = 0.;

   for ( int i=1;i<101;i++){
     Double_t b=blow + i*bstep;
     sum1 = sum1+b;
     sum2 = sum2 + b*(Npart(Ap,At,b)+Npart(At,Ap,b));  // Gosset function
   }

   Double_t NGosset = sum2/sum1;
   Double_t xa = Ene;    //  use TAPS data now (PRC 56 (1997) R2920)
   
   Double_t sigpi0;
   Double_t sigother;
   Double_t sigeta;
   Double_t ppi0;
   Double_t peta;
   Double_t pother;
   Double_t pom;

   if ((Part==7)||(Part==8)||(Part==9)){
     Double_t spi0CC = exp(5.659+1.706*log(xa)-0.620*pow(log(xa),2)+
                     0.0511*pow(log(xa),3));
                                                                        // pi0

     Double_t spi0CaCa = exp(7.429+1.738*log(xa)-.5599*pow(log(xa),2)-
                       .00210*pow(log(xa),3));
     Double_t alfapi0=0.75+0.045*xa; // linear alpha used instead of constant


     Double_t kpi0 =1.05*sqrt(spi0CC*spi0CaCa)/pow(12*40,alfapi0);

     Double_t spi0ApAt = kpi0*pow((Ap*At),alfapi0); 
                                       //extrapolate pi0 cross sec. (in mb)
     ppi0 = spi0ApAt/(sigmax*Apart);
     sigpi0 = ppi0*NGosset*sigreac;    // in mb
     if (Part==7) {
       pom = sigpi0;
     } else if (Part==8) {
       pom = sigpi0*pip; 
     } else {
       pom = sigpi0*pim; 
     }

   } else {

     Double_t setaCC = exp(0.470+4.860*log(xa)-1.321*pow(log(xa),2)); 
                                                                // eta C + C
     Double_t setaCaCa = exp(2.6885+4.2342*log(xa)-1.774*pow(log(xa),2)); 
                                                                // eta Ca + Ca
     Double_t alfaeta = 1.0757 - 0.19524*xa;

     Double_t keta = 1.05*sqrt(setaCC*setaCaCa)/pow(12*40,alfaeta); 
     Double_t setaApAt = keta*pow((Ap*At),alfaeta);   
                                      //  extrapolate eta cross sec. (in mb)
     peta = setaApAt/(sigmax*Apart);
     sigeta = peta*NGosset*sigreac;   // in mb
  
     pother = peta*ratiosignew(T,(1000*makeStaticData()->GetParticleMass(Part)));  
                                      // mt-scale eta for other particles
     sigother = pother*NGosset*sigreac; 
     if (Part==17) {
       pom = sigeta;
     } else {
       pom = sigother;
     }  
   }
   return pom;  // return cross section in mb
}


Double_t PCross::cross(char * Part, Double_t blow, Double_t bup) {
  return cross(makeStaticData()->GetParticleID(Part),blow,bup);
}


Double_t PCross::cross(Int_t Part, Double_t blow, Double_t bup) {
  //interface to function calc(*x,*par) called by TF1
  Double_t x[1];
  Double_t par[7];
  if (check(Part)) {
    x[0] = Ebeam;
    par[0] = Part;
    par[1] = blow;
    par[2] = bup;
    par[3] = AT;
    par[4] = ZT;
    par[5] = AP;
    par[6] = ZP;
    return calc(x,par);
  }
  return 0.;
}


Bool_t PCross::check(Int_t Part) {  // check validity of particle number
  if (!sys) {
    cout << "Sorry, first you have to setSystem()!"  << endl;
    return kFALSE;
  }

  if (Part==0) {
    cout << "List of presently supported particles:"<< endl;
    cout << "7)   pi0"<< endl;
    cout << "8)   pi+"<< endl;
    cout << "9)   pi-"<< endl;
    cout << "10)  K0L"<< endl;
    cout << "16)  K0S"<< endl;
    cout << "17)  eta"<< endl;
    cout << "52)  omega"<< endl;
    cout << "55)  phi"<< endl;
    return kFALSE;

  } else if ((Part==7)||(Part==8)||(Part==9)||(Part==10)||(Part==16)||
             (Part==17)||(Part==52)||(Part==55)) {
    return kTRUE;

  } else {
    cout <<"Sorry, this particle is not yet supported!" << endl;
    return kFALSE;
  }
}

Double_t PCross::calcT(Double_t Ebeam) { // fireball temperature
//  T = 23.+47.*Ebeam -7.*Ebeam*Ebeam;   // in MeV (old TAPS formula)
//  chemical freeze-out according to Ref 3 (Cleymans et al. PRC 73 (2006) 034905)
    Double_t sqs = sqrt(2.*0.939*(Ebeam + 2.*0.939));
    Double_t mu = 1.308/(1.+0.273*sqs);  // chemical potential in GeV
    return 0.166 - 0.139*mu*mu - 0.053*mu*mu*mu*mu;   // T in GeV
}

