/**************************************************
**      This class is licensed under LGPL        **
**  - free for personal and commercional use -   **
**                                               **
**   ---------------------------------------     **
**   Defining the fermi motion distributions     **
**   ---------------------------------------     **
**                                               **
**   Author: L. Witthauer & M. Dieterle 2009     **
**                  Version 1.0                  **
**                                               **
**************************************************/

#include "PFermiDistributions.h"

//------------------------------------------------------------------------

PFermiDistributions::PFermiDistributions() {}

PFermiDistributions::PFermiDistributions(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de, key)  {

    DeutC = new double[13];
    DeutD = new double[13]; 
    DeutM2 = new double[13];
   
    //Now read my own target:
    //cout << is_pid << endl;
 
    TString Target(makeStaticData()->GetParticleName(is_pid));
    Double_t target_mass = makeStaticData()->GetParticleMass(is_pid);

#if 1
    Target.ToUpper();
   

    if (Target == "D")
    {
       SetDeuteronValues();
       version_flag = NUCLEAR_FERMI_D;
       //       FM = new TF1("FM",PFermiDistributions::FermiDisD, 0.,1.,2);
       MeV2GeV = 1.;
    }
    else if (Target == "HE3")
    {
	version_flag = NUCLEAR_FERMI_3HE;
	//FM = new TF1("FM",PFermiDistributions::FermiDis3He, 0.,1000.,2);
	MeV2GeV = 0.001;
    }   
    else if (Target == "ALPHA")
    {
	version_flag = NUCLEAR_FERMI_4HE;
	//FM = new TF1("FM",PFermiDistributions::FermiDis4He, 0.,1.,1);
       MeV2GeV = 1.;
    }
    else if (Target == "7LI")
    {
	version_flag = NUCLEAR_FERMI_7LI;
	//FM = new TF1("FM",PFermiDistributions::FermiDis7Li, 0.,1000.,3);
       MeV2GeV = 0.001;
    }
    else if (Target == "12C")
    {
	version_flag = NUCLEAR_FERMI_12C;
        //FM = new TF1("FM",PFermiDistributions::FermiDis12C,0.,1000.,2);
       MeV2GeV = 0.001;
    }
    else if (Target == "40CA" || Target == "93NB" || Target == "208PB")
    {
	version_flag = NUCLEAR_FERMI_40CA;
	//FM = new TF1("FM",PFermiDistributions::FermiDis40Ca,0.,1000.,2);
       MeV2GeV = 0.001;
    }  
    else if (target_mass > 0.5*(makeStaticData()->GetParticleMass("40Ca")
				+makeStaticData()->GetParticleMass("12C"))) {
	Warning("SetTargetMaterial","Target Material %s not found, assuming 40Ca",
		makeStaticData()->GetParticleName(is_pid)); 
	version_flag = NUCLEAR_FERMI_40CA;
	MeV2GeV = 0.001;
    } 
    else if (target_mass > 0.5*(makeStaticData()->GetParticleMass("7Li")
				+makeStaticData()->GetParticleMass("12C"))) {
	Warning("SetTargetMaterial","Target Material %s not found, assuming 12C",
		makeStaticData()->GetParticleName(is_pid)); 
	version_flag = NUCLEAR_FERMI_12C;
	MeV2GeV = 0.001;
    } 
    else if (target_mass > 0.5*(makeStaticData()->GetParticleMass("7Li")
				+makeStaticData()->GetParticleMass("alpha"))) {
	Warning("SetTargetMaterial","Target Material %s not found, assuming 7Li",
		makeStaticData()->GetParticleName(is_pid)); 
	version_flag = NUCLEAR_FERMI_7LI;
	MeV2GeV = 0.001;
    } 
    else {
	version_flag = 0;
	Error("SetTargetMaterial","Target Material not found"); 
    }
#endif
};


PDistribution* PFermiDistributions::Clone(const char*delme) const {
    return new PFermiDistributions((const PFermiDistributions &)* this);
};

PFermiDistributions::~PFermiDistributions(){

    //  delete FM;
   delete [] DeutC;
   delete [] DeutD;
   delete [] DeutM2;

};

Double_t  PFermiDistributions:: EvalPar(const double* x, const double* par) {
    return Eval(x[0]);
}

Double_t PFermiDistributions::Eval(Double_t x, Double_t y , Double_t z , Double_t t ) const {
    // Returns connected fermi distribution
    Double_t p=x/MeV2GeV;
    if (version_flag == NUCLEAR_FERMI_D)
	return FermiDisD(&p,NULL);
    else if (version_flag == NUCLEAR_FERMI_3HE)
	return FermiDis3He(&p,NULL);
    else if (version_flag == NUCLEAR_FERMI_4HE)
	return FermiDis4He(&p,NULL);
    else if (version_flag == NUCLEAR_FERMI_7LI)
	return FermiDis7Li(&p,NULL);
    else if (version_flag == NUCLEAR_FERMI_12C)
	return FermiDis12C(&p,NULL);
    else if (version_flag == NUCLEAR_FERMI_40CA)
	return FermiDis40Ca(&p,NULL);
    return 0.;
}


//---------------------- Fermi Distribution of Deuteron  -----------------------
double PFermiDistributions::FermiDisD( double *x, double *par ) const {

   const double pi    = 3.1415927;
   const double hbar  = 0.197463569747999998;  // unit: GeV*fm/c 
   double par0 = 0.;
   double par1 = 0.;
   for (int j = 0; j < 13; j++)
   {
      par0 += DeutC[j]/(x[0]*x[0]/(hbar*hbar) + DeutM2[j]);
      par1 += DeutD[j]/(x[0]*x[0]/(hbar*hbar) + DeutM2[j]);
   }
   par0  = TMath::Sqrt(2./pi)*par0;
   par1  = TMath::Sqrt(2./pi)*par1;
   Double_t mom = (par0*par0+par1*par1)*x[0]*x[0]/(hbar*hbar);
  
   return mom;   
}
//---------------------- Fermi Distribution of Helium-3 ---------------------

double PFermiDistributions::FermiDis3He( double *x, double *par ) const {

  const double pi = 3.1415927;
  const double a  = 7.09078;
  const double b  = 5.38753;
  const double c  = 9.90202;
  const double d  = 0.779408;

  Double_t par0 = (4./TMath::Sqrt(pi))*TMath::Power(a,3./2.);
  Double_t par1 = x[0]*x[0]*25./(1.E6);
  Double_t mom = d*par0*par1*(TMath::Exp(-par1*a) + c*TMath::Exp(-TMath::Sqrt(par1)*b));

  return mom;
}
//---------------------- Fermi Distribution of Helium-4 ---------------------

double PFermiDistributions::FermiDis4He( double *x, double *par )  const {

   const double hbar  = 0.197463569747999998;
   const double a = 0.7352;
   const double b = 0.05511;

   Double_t par0 = TMath::Exp(-x[0]*x[0]/(hbar*hbar*a));
   Double_t mom = x[0]*x[0]/(hbar*hbar*b)*par0;

   return mom;

}
//---------------------- Fermi Distribution of Lithium-7 ---------------------
double PFermiDistributions::FermiDis7Li( double *x, double *par ) const {

   const double a = 1.2e-4;
   const double b = 6.87e-3;
   const double c = 110.;
   
   Double_t par0 = x[0]*x[0]*a;
   Double_t par1 = (1. + b*x[0]);
   Double_t par2 = x[0]/c;
   Double_t mom  = par0*par1*TMath::Exp(-par2*par2);

   return mom;

}
//---------------------- Fermi Distribution of Carbon-12 ---------------------
double PFermiDistributions::FermiDis12C( double *x, double *par ) const {

   const double pi    = 3.1415927;
   const double a = 1/0.416;
   const double b = 1/0.23;
   const double c = 0.04;
   
   Double_t par0 = (4./TMath::Sqrt(pi))*TMath::Power(a,3./2.);
   Double_t par1 = x[0]*x[0]*25./(1.E6);
   Double_t mom  = par0*par1*(TMath::Exp(-par1*a) + c*TMath::Exp(-TMath::Sqrt(par1)*b));

   return mom;

}

//---- Fermi Distribution of Calcium-40 (can be used for heavier nuclei) -----------------
double PFermiDistributions::FermiDis40Ca( double *x, double *par ) const {

   const double pi    = 3.1415927;
   const double a = 1/0.42;
   const double b = 1/0.23;
   const double c = 0.04;
   
   Double_t par0 = (4./TMath::Sqrt(pi))*TMath::Power(a,3./2.);
   Double_t par1 = x[0]*x[0]*25./(1.E6);
   Double_t mom  = par0*par1*(TMath::Exp(-par1*a) + c*TMath::Exp(-TMath::Sqrt(par1)*b));

   return mom;

}

//------ Calculate Values for Deuteron Wave Function ----------------
void PFermiDistributions::SetDeuteronValues( ){
   
   const double alpha = 0.23162461;            // unit: 1/fm
   const double mZero = 1.;                    // unit: 1/fm
   double cj[13]={  0.88688076E0,              // unit: 1/sqrt(fm)
		   -0.34717093E0,
		   -0.30502380E1,
		    0.56207766E2,
		   -0.74957334E3,
		    0.53365279E4,
		   -0.22706863E5,
		    0.60434469E5,
		   -0.10292058E6,
		    0.11223357E6,
		   -0.75925226E5,
		    0.29059715E5,
		    0.0};
   double dj[13]={  0.23135193E-1,             // unit: 1/sqrt(fm)
		   -0.85604572E0,
		    0.56068193E1,
		   -0.69462922E2,
		    0.41631118E3,
		   -0.12546621E4,
		    0.12387830E4,
		    0.33739172E4,
		   -0.13041151E5,
		    0.19512524E5,
		    0.0,
		    0.0,
		    0.0};
   double sum1 = 0.;
   double sum2 = 0.;
   double sum3 = 0.;
   double rtmp;
   double mj[13];
   double mj2[13];
   int n;
   int n1;
   int n2;
   int temp;

   for ( int j = 0; j < 13; j++)
   {
      mj[j]  = alpha + j*mZero;
      mj2[j] = mj[j]*mj[j];
   } 
   
   cj[12] = 0.;
   for ( int i = 0; i < 12; i++) cj[12] = cj[12] - cj[i];
   
   for ( int k = 0; k < 5; k++)
   {
      rtmp = dj[k*2]/mj2[k*2] + dj[k*2+1]/mj2[k*2+1];
      sum1 = sum1 + rtmp;
      rtmp = dj[k*2] + dj[k*2+1];
      sum2 = sum2 + rtmp;
      rtmp = dj[k*2]*mj2[k*2] + dj[k*2+1]*mj2[k*2+1];
      sum3 = sum3 + rtmp;
   }

   n  = 12;
   n1 = 11;
   n2 = 10;

   for ( int l = 0; l < 3 ; l++)
   {  
      dj[n2] = -mj2[n1]*mj2[n]*sum1 + (mj2[n1]+mj2[n])*sum2 - sum3;
      dj[n2] = dj[n2] * mj2[n2]/(mj2[n]-mj2[n2])/(mj2[n1]-mj2[n2]); 

      temp = n2;
      n2 = n1;
      n1 = n;
      n = temp; 
   }

   for (int j=0;j<13;j++)
   {
      DeutC[j]   = cj[j];
      DeutD[j]   = dj[j];
      DeutM2[j] = mj2[j];
   }
}


//-----------------------------------------------------------------------
ClassImp(PFermiDistributions);

