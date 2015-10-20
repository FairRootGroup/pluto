////////////////////////////////////////////////////////
//  Momentum smearer implementation file
//
//  
//
////////////////////////////////////////////////////////


#include "PHadesParticleSmearer.h"
#include <iostream>
#include <cmath>

using namespace std;

 
PHadesParticleSmearer:: PHadesParticleSmearer() {
    resolution_factor=1;
    fPriority = FILTER_PRIORITY+1; //for HADES it is standard to filter first and then smear
}

bool PHadesParticleSmearer::Modify(PParticle ** stack, int *decay_done, int * num, int stacksize) {


    for (int i=0; i< *num; i++) {
	PParticle * cur = stack[i];
//	cur->Print();
	

	if (cur->Is("pi+") || cur->Is("pi-")|| cur->Is("p") 
	    || cur->LeptonN()) {

/*
		Double_t Impulsx = cur->Px();
		//Impulsx = 5.*Impulsx;

*/
		//double mom_measured = PUtils::sampleGaus(Impulsx, Impulsx*0.1);
		//cur->SetPx(mom_measured);

		double res = GetResolution(cur,4); 
       		Smear(cur, res);// momentum smearing / 4mdcs / geidar's parametrisation 


 	    cur->ResetE();

	}
	
    }
    
    
    return kTRUE;
};




double PHadesParticleSmearer::GetResolution(PParticle *pPart, Int_t iSetup)
{
   //double resolution_factor = 3;
   double theta = TMath::RadToDeg()*(pPart->Theta());  // polar angle in degree
   double phi = 90.-TMath::RadToDeg()*(pPart->Phi());  // azimuthal angle in degree
   double mass = pPart->M()*1000.;
   double p = (pPart->Vect()).Mag()*1000.;
   if (phi<0.) phi += 360.;                // TLorentzVector uses ATan2(x,y)
   if(theta <= 15 || theta > 85) return 0.004;
   Int_t thetaBin = int((theta-15)/10);
   Int_t phiBin = int((phi - int(phi/60)*60)/12);
   double c1, c2, res;
   switch(iSetup)
      {
      case 2:
         c1 = detkick[thetaBin][phiBin];
         res = c1*p*p;
         break;

      case 3:
         c1 = detrk3[thetaBin][phiBin];
         c2 = multscrk3[thetaBin][phiBin];
         res = sqrt(c1*c1*p*p*p*p + c2*c2*(p*p + mass*mass));
         break;

      case 4:
         c1 = detrk4[thetaBin][phiBin];
         c2 = multscrk4[thetaBin][phiBin];
         res = sqrt(c1*c1*p*p*p*p + c2*c2*(p*p + mass*mass));
         break;

      default:
         res = 100;
      }
      res=resolution_factor*res/1000;
      return res;
}



void PHadesParticleSmearer::Smear(PParticle *p, double gamma) 
{
  //double ptot=(p->Vect()).Mag();
  double ptot=p->P();
  if (ptot<=0.) return;
  double th=p->Theta(), ph=p->Phi(), sth =sin(th), mass=p->M();
  if (fabs(sth)<1.e-6) sth = (sth>0.) ? 1.e-6 : -1.e-6;
  double smeared_ptot = ( rand.Gaus(ptot, gamma ) );
  if(gamma==0.1) smeared_ptot = ptot+0.04;
//////
//  smeared_ptot= ptot;
//cout<<ptot<<"  "<<gamma<<endl;
//////
//if(fabs(ptot-smeared_ptot)>0.15) cout << ptot << " " << smeared_ptot << endl;
  double smeared_theta = rand.Gaus(th, 0.0023 );
  double smeared_phi = ( rand.Gaus(sth*ph, 0.0023 ) )/sth;
  double theta_check=smeared_theta-2.*TMath::Pi()*int(smeared_theta/(2.*TMath::Pi()));
  if (theta_check>TMath::Pi()||theta_check<0.) smeared_phi -= TMath::Pi();
  p->SetRho(smeared_ptot);                    // keeping theta, phi constant
  p->SetTheta(smeared_theta);                 // keeping |p|, phi constant
  p->SetPhi(smeared_phi);                     // keeping |p|, theta constant
  double newp=p->P();
  p->SetE(sqrt(newp*newp+mass*mass));         // reset energy
}

// Geidar parametrisation
double PHadesParticleSmearer::detrk4[7][5]     = {{3.04669e-06, 2.92245e-06, 2.88894e-06, 2.54182e-06, 1.79336e-06},
                                       {3.45534e-06, 3.52212e-06, 3.43443e-06, 3.3966e-06,  3.31397e-06},
                                       {5.01234e-06, 5.46206e-06, 5.25006e-06, 5.06055e-06, 4.87875e-06},
                                       {6.29485e-06, 7.83364e-06, 8.64567e-06, 7.44434e-06, 6.44697e-06},
                                       {5.90814e-06, 7.49372e-06, 8.44064e-06, 7.33158e-06, 5.8983e-06} ,
                                       {5.00201e-06, 6.82347e-06, 8.42634e-06, 6.64071e-06, 4.32283e-06},
                                       {4.57113e-06, 7.85962e-06, 9.27056e-06, 7.13838e-06, 4.2649e-06 }};

double PHadesParticleSmearer::multscrk4[7][5] = {{0.00712136, 0.00774792, 0.00745848, 0.0072562, 0.0077318},
                                      {0.00947987, 0.00919469, 0.00923363, 0.0090577, 0.0088233},
                                      {0.0104807 , 0.0105168 , 0.0114445 , 0.0109526, 0.0103268},
                                      {0.0122493 , 0.013979  , 0.0148295 , 0.0142261, 0.0125651},
                                      {0.0132593 , 0.0151324 , 0.0171913 , 0.0157784, 0.0133762},
                                      {0.0127184 , 0.0166672 , 0.0197425 , 0.0166028, 0.0127183},
                                      {0.0121717 , 0.0187349 , 0.0235083 , 0.0183666, 0.0114857}};

double PHadesParticleSmearer::detrk3[7][5] = {{5.08616e-06, 5.57472e-06, 5.90756e-06, 5.30146e-06, 5.12577e-06},
                                   {8.64475e-06, 8.44705e-06, 8.58983e-06, 7.7759e-06 , 6.88931e-06},
                                   {1.23261e-05, 1.42034e-05, 1.55287e-05, 1.31425e-05, 1.10076e-05},
                                   {1.75029e-05, 2.59217e-05, 2.72391e-05, 2.46139e-05, 1.84151e-05},
                                   {1.75554e-05, 2.25547e-05, 2.5371e-05 , 2.30453e-05, 1.8036e-05 },
                                   {1.5159e-05 , 2.14014e-05, 2.37335e-05, 2.17028e-05, 1.56393e-05},
                                   {1.70025e-05, 2.42703e-05, 2.65468e-05, 2.33554e-05, 1.5726e-05 } };

double PHadesParticleSmearer::multscrk3[7][5] = {{0.00705671, 0.00707791, 0.00719795, 0.00686656, 0.00468438},
                                     { 0.00791494, 0.00971418, 0.00922887, 0.00921739, 0.00882323},
                                      {0.0101972 , 0.0108932 , 0.0109951 , 0.0111292 , 0.0103408 },
                                     { 0.0136591 , 0.0147924 , 0.0152897 , 0.0153152 , 0.0135126 },
                                      {0.014621  , 0.0166669 , 0.0175123 , 0.0176646 , 0.0157295 },
                                      {0.0146804 , 0.0180094 , 0.0183208 , 0.0172132 , 0.014879 } ,
                                      {0.0130367 , 0.018929  , 0.0221153 , 0.0189616 , 0.0147333 }   };

double PHadesParticleSmearer::detkick[7][5] = {{0.000095203, 0.000095204, 0.000104141, 0.000099331, 0.000099331},
                                    {0.000151284, 0.000148333, 0.000151814, 0.000147953, 0.000143505},
                                    {0.000197612, 0.000196251, 0.000191317, 0.000195112, 0.000203217},
                                    {0.000156319, 0.000157184, 0.000163081, 0.000157349, 0.000162052},
                                    {0.00014439 , 0.000153818, 0.000170083, 0.000160573, 0.000139122},
                                    {0.000199936, 0.000236463, 0.000262118, 0.000240549, 0.000199932},
                                    {0.000226558, 0.000296774, 0.00032661 , 0.000292281, 0.000226592}};






ClassImp(PHadesParticleSmearer)
