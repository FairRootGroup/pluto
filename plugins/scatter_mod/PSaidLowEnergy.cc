///////////////////////////////////////////////////////
//    Pluto SAID parametrization implementation file
//    Low energy SAID pp-elastic parametrization.
//    Created by Nodar Lomidze and Mirian Tabidze;  
//    ANKE/PAX Collaboration 03.09.09  
//
//    For each theta (0-180) deg. dsg/dT_kin distributions was fitted using function
//    a*pow[T_kin,b*pow(T_kin,c)] + d, where a,b,c, and d parameters are functions of 
//    theta, T_kin is proton kinetic energy.
//    For approximation this parameters full theta region was divided into 5 intervals 
//    (0-1), (1-4), (4-20), (20-40) and (40-90) deg.
//    In each interval polynomial functions (pol4-pol6) or function 
//    pol1+A*x^B*sin^C(x)*exp(D*x)  was used (see code below), 
//    where A,B,C and D are parameters.
//    In (90-180) deg. interval same function was used as in (90-0) deg. 
//
//    Differences from SAID above 4 MeV incident proton kinetic energy in theta
//    (25-90) deg. are less then 1%, in (4-25) deg. ~4-6% and in (1-4) deg. (10-15)%. 
//    Below 4 MeV parametrization is worse ~(10-20)% in all theta intervals.
//
//
////////////////////////////////////////////////////////


#include "PSaidLowEnergy.h"

const long double d2r=3.14159265358979323846/180.;



PSaidLowEnergy::PSaidLowEnergy() {
} ;

PSaidLowEnergy::PSaidLowEnergy(const Char_t *id, const Char_t *de) :
    PAngularDistribution(id, de) {
    dim=90-1;

    y.Set(dim+1);
    aa.Set(dim);
    told=-1.;

} ;

PDistribution* PSaidLowEnergy::Clone(const char*delme) const {
    return new PSaidLowEnergy((const PSaidLowEnergy &)* this);
};


Bool_t PSaidLowEnergy::IsValid(void) {
    
    if (!direct_sampling_done) {
	Fatal("IsValid","Sampling not finished");
    }
    return kTRUE;
};


Bool_t PSaidLowEnergy:: Init(void) {
    if (!PAngularDistribution:: Init()) return kFALSE;


    if (!beam || !target) {
	Warning("Init","beam or target not found");
	return kFALSE;
    }

    return kTRUE;
}


Bool_t PSaidLowEnergy::Prepare(void) {
    return kTRUE;
}

double PSaidLowEnergy::SamplePolarAngle(double r) {
    //boost particle 2 such to be in particle 1 rest frame

    PParticle tmp(beam);  //particle under investigation. Make better a copy
    tmp.Boost(-target->BoostVector());
    Double_t tlab = tmp.KE();

    return sample(r,tlab);
}

//mita begin
/*
#include <TCanvas.h>
#include <TH1F.h>
int I_print=0;
TCanvas *cc = new TCanvas("c","SAID parametrisation",200,10,800,700);
TH1F *hf0t=new TH1F("hf0t",";#Theta_{lab} (deg.)",180,0,90);
*/
//mita end


double PSaidLowEnergy::dsdw(double th, double tlab) {
    // normalized pp elastic differential cross sections (mb/sr)
    // arguments: scattering angle (deg), lab beam kinetic energy (GeV)

    double Tlab = tlab*1000;    

    double a0=825.770315997777;
    double a1=-36.1806390269961;
    double a2=0.736508337604257;
    double a3=-0.00245475650713848;
    double a4=-9.46167840068291e-05;
    double a5=1.18647782825118e-06;
    double a6=-4.23920644982741e-09;
    
    double a01=158509949.824459;
    double a11=-1.05462493339422;
    double a21=0.999464445523888;
    double a31=-0.331696205972712;
    double a41=-1.9817093971156;
    double a51=339.459572120675;
    
    double a02=7013162.6023937;
    double a12=-2.60581893459221;
    double a22=-1.17622988844922 ;
    double a32=-0.0572294874172003;
    double a42=-242.76957186148;
    double a52=4972.40753398661;
    
    double a03=7363793.10073931;
    double a13=-2.78666712844778;
    double a23=-1.19175331577741;
    //    double a33=-0.0168240994020085;
    double a43=2919.73115823142;
    double a53=-11247.1848613229;
    
    double b0=-1.11987335717749;
    double b1=0.0375681831247435;
    double b2=-0.000806720177542215;
    double b3=7.31226869503504e-06;
    double b4=-2.40772264489991e-08 ;
    
    double b01=5.49631512561751;
    double b11=-1.41153243358003;
    double b21=0.0873325395997571;
    double b31=-0.00210915057288232;
    double b41=1.78555710528402e-05;
    
    double b02=-2.01651356175559;
    double b12=0.00952134681173633;
    double b22=-0.00424259230244138;
    double b32=0.000815075673344715;
    double b42=-6.09637287432679e-05;
    double b52=1.5439817446594e-06;
    
    double c0=-0.0899947020336036;
    double c1=0.0165457502998486;
    double c2=-0.00035225623300735;
    double c3=3.18377219079531e-06;
    double c4=-1.04867548784599e-08;
    
    double c01=4.02283676522501;
    double c11=-0.571280456418181;
    double c21=0.0286964703485152;
    double c31=-0.000606564121973199;
    double c41=4.6543463351715e-06;
    
    double c02=0.000449138526079407;
    double c12=-0.00323470722989942;
    double c22=0.00101229269064745;
    double c32=-8.04320221005051e-05;
    double c42=1.38972236339651e-06;
    double c52=1.51442646510768e-08;
    
    double d0=0.14637558766108;
    double d1=0.242429384264589;
    double d2=-0.00515619107849952;
    double d3=4.5272537785055e-05;
    double d4=-1.42514150093272e-07;
    
    double d01=90.5702938678518;
    double d11=-12.623103132398;
    double d21=0.63864747288047;
    double d31=-0.0137296412910617;
    double d41=0.000107543327455007;
    
    double d03=95.4776480269654;
    double d13=-2.46798492423323;
    double d23=-1.27436869909584;
    double d33=-0.478148940786103;
    double d43=11.6930721555554;
    double d53=-64.5772975461858;
    double d63=-0.702778910136937;
    double d73=0.014105237506306;
    
    double as, bs, cs,ds;
    
    if(th>40) {
      as=a0+a1*th+a2*pow(th,2)+a3*pow(th,3)+a4*pow(th,4)+a5*pow(th,5)+a6*pow(th,6);
      bs=b0+b1*th+b2*pow(th,2)+b3*pow(th,3)+b4*pow(th,4);
      cs=c0+c1*th+c2*pow(th,2)+c3*pow(th,3)+c4*pow(th,4);
      ds=d0+d1*th+d2*pow(th,2)+d3*pow(th,3)+d4*pow(th,4); 
    }

    if(th>20 && th <= 40) {
      as=a01*pow(th,a11)*pow(sin(th*d2r),a21)*exp(a31*th)+a41*th+a51;
      bs=b01+b11*th+b21*pow(th,2)+b31*pow(th,3)+b41*pow(th,4);
      cs=c01+c11*th+c21*pow(th,2)+c31*pow(th,3)+c41*pow(th,4);
      ds=d01+d11*th+d21*pow(th,2)+d31*pow(th,3)+d41*pow(th,4);  
    }
    
    if(th> 3.99 && th<=20) {
      as=a02*pow(th,a12)*pow(sin(th*d2r),a22)*exp(a32*th)+a42*th+a52;
      bs=b02+b12*th+b22*pow(th,2)+b32*pow(th,3)+b42*pow(th,4)+b52*pow(th,5);
      cs=c02+c12*th+c22*pow(th,2)+c32*pow(th,3)+c42*pow(th,4)+c52*pow(th,5);
      ds=d03*pow(th,d13)*pow(sin(th*d2r),d23)*exp(d33*th)+d43*th+d53+d63*pow(th,2)+d73*pow(th,3);
    }
 
    if(th> 1 && th<=3.99) {
      as=a03*pow(th,a13)*pow(sin(th*d2r),a23)*exp(a32*th)+a43*th+a53;
      bs=b02+b12*th+b22*pow(th,2)+b32*pow(th,3)+b42*pow(th,4)+b52*pow(th,5);
      cs=c02+c12*th+c22*pow(th,2)+c32*pow(th,3)+c42*pow(th,4)+c52*pow(th,5);
      ds=d03*pow(th,d13)*pow(sin(th*d2r),d23)*exp(d33*th)+d43*th+d53+d63*pow(th,2)+d73*pow(th,3);
    }
    
    if(th<=1) {
      as=(901707552-762)*th+762;
      bs=(-2.01+0.628)*th-0.628;
      cs=(-0.0024+0.229)*th-0.229;
      ds=(10244.8+267.38)*th-267.38;
    }
    
    double sinth = sin(th*d2r);

 /*
      double dSdW, dSdT;
      I_print = I_print + 1;
      dSdW  = as*pow(Tlab,bs*pow(Tlab,cs)) + ds;
      dSdT  = sinth*dSdW;
      if(I_print <= 110) {
      cout << I_print << " Tkin = " << Tlab << endl;
      cout << " Theta " << th << " dS/dW " << dSdW 
      << " dS/dT " << dSdT << endl;
      }
      hf0t->Fill(th);
*/
    
    return sinth*(as*pow(Tlab,bs*pow(Tlab,cs)) + ds);

}

double PSaidLowEnergy::sample(double r, double tlab) {
    const double fac=1.3;
    int i, j;
    // _____________________________________________________________________
    // enter whenever beam energy changes, to calculate integrals and slopes
    if (told!=tlab) {
	told=tlab;
	double sum=0.;
	for (i=0;i<=dim;++i) y[i]=fac*dsdw((double)(i+1),tlab);
	for (i=0;i<dim;++i) {
	    sum += (y[i+1]>y[i]) ? y[i+1] : y[i] ;
	    aa[i]=sum;
	}
    }
  
    // _____________________________________________________________________
    // sample cm scattering angle with the rejection method
    double rr, f, tf, r2, c0, r1=r;
 repeat:
    i=(r1>.5);
    rr=(2.*r1-i)*aa[dim-1];
    for (j=0;j<dim;++j) {if (rr<aa[j]) break;}
    if (j) rr-=aa[j-1];
    tf = (y[j+1]>y[j]) ? y[j+1]: y[j];
    c0=rr/tf+1.+j;
    f=dsdw(c0,tlab);
    //r2=REngine->Rndm();
    r2=PUtils::sampleFlat();
    if (r2>f/tf) goto fail;
    if (f<=tf) return (1-2*i)*cos(c0*d2r);

    Warning("sample","sampling failed f=%lf, tf=%lf",f,tf);
    exit(1);

 fail:
    //r1=REngine->Rndm();
    r1=PUtils::sampleFlat();
    goto repeat;
}


ClassImp(PSaidLowEnergy)
