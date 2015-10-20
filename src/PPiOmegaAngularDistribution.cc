/////////////////////////////////////////////////////////////////////
//
//  Angular distributions in pi+p -> omega + X
//
//                                  Author: Johan Messchendorp
//                                  Reimplemented I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PPiOmegaAngularDistribution.h"

const long double pi=3.14159265358979323846;

ClassImp(PPiOmegaAngularDistribution)

PPiOmegaAngularDistribution::PPiOmegaAngularDistribution() {

} ;

PPiOmegaAngularDistribution::PPiOmegaAngularDistribution(const Char_t *id, const Char_t *de) :
    PAngularDistribution(id, de) {
    beam = NULL;
    target = NULL;
    e_cm=0;
    version=1;
} ;

PDistribution* PPiOmegaAngularDistribution::Clone(const char*delme) const {
    return new PPiOmegaAngularDistribution((const PPiOmegaAngularDistribution &)* this);
};




Bool_t PPiOmegaAngularDistribution::IsValid(void) {
    
    if (!direct_sampling_done) {
	Fatal("IsValid","Sampling not finished");
    }
    
    return kTRUE;
};

Bool_t PPiOmegaAngularDistribution:: Init(void) {
    if (!PAngularDistribution:: Init()) return kFALSE;

    //In addition we need the beam and target
    //done already in PAngularDistribution

    if (!beam || !target) {
	Warning("Init","beam or target not found");
	return kFALSE;
    }

    return kTRUE;
}

Bool_t PPiOmegaAngularDistribution:: Prepare(void) {
    getPiN_wN_param();
    return kTRUE;
}


void PPiOmegaAngularDistribution:: getPiN_wN_param() {
  // energy-dependent parameters of the cm angular distribution 
  // for w production in pi+N --> N+w
  if (parent->M()!=e_cm) e_cm=parent->M();
  else return;                         // no change in values

  const double y[]={exp(-8.),exp(-4.),exp(-3.2),exp(-2.4),1.},
    dx[]={1.,.2,.2,.6}, dy[]={y[1]-y[0],y[2]-y[1],y[3]-y[2],y[4]-y[3]};
  int i;
  double h=.0069*exp(2.873*e_cm), a=0.;   // h parametrized by J.Messchendorp
  PiN_w_y[0]=1.+h*y[0];
  PiN_w_h=h;
  for (i=0;i<4;++i) {
    PiN_w_y[i+1]=1.+h*y[i+1];
    PiN_w_slope[i]=h*dy[i]/dx[i];
    a+=dx[i]*(0.5*h*dy[i]+PiN_w_y[i]);
    PiN_w_area[i]=a;
  }
}


double PPiOmegaAngularDistribution::SamplePolarAngle(double r) {
  const double xj[5]={-1.,0.,.2,.4,1.};
  double rr, f=-1, tf=-1, r2, c0=-1, r1=r, a=-1, b=-1, area=-1, x1, x2, delta;
  int i, fl=0 ;

 again:
  if (version==PI_OMEGA_piNNw) {                   // pi + N --> N + w, PRD2 (1970) 2564
    fl=2;
    rr=r1*PiN_w_area[3];
    for (i=0;i<4;++i) if (rr<PiN_w_area[i]) break;
    area=(i)?PiN_w_area[i-1]:0.;
    x1=PiN_w_y[i];
    x2=x1/PiN_w_slope[i];
    delta=1.-2.*(area-rr)/(x1*x2);
    c0=-x2*(1.-sqrt(delta))+xj[i];
    f=1.+PiN_w_h*exp(-4.+4.*c0);
    tf=(c0-xj[i])*PiN_w_slope[i]+x1;
  } else if (version==PI_OMEGA_piPPpiw) {                 // pi+ + P --> pi+ P + w,
                                        // Alff PR 145 (1966) 361
    fl=3;
    c0 = log((r1*21.1545+3.31)/9.);
    f  = 3. + exp(3.*c0);
    tf = 9. * exp(c0); c0 = -c0;
  } else if (version==PI_OMEGA_piPDw) {                   // pi+ + P --> pi+ P + w,
    fl=4;                               // Alff PR 145 (1966) 361
    c0 = log((r1*23.504+3.679)/10.);
    f  = 2. + exp(2.3*c0);
    tf = 10. * exp(c0); c0 = -c0;
  }

  else return c0=PUtils::sampleFlat()*2.-1.; // leave it isotropic for the moment
    
  r2=PUtils::sampleFlat();
  if (r2>f/tf) {                        // reject
    r1=PUtils::sampleFlat();
    goto again;
  }
  if (tf<f) {                           // diagnostics
    printf("PPiOmegaAngularDistribution::samplePolarAngle: algorithm failed %d\n",fl);
    printf("a=%f b=%f: c0=%f tf=%f f=%f\n",a,b,c0,tf,f);
    r1=PUtils::sampleFlat();
    goto again;
  }
  return c0;                            // accept
    


}

void PPiOmegaAngularDistribution::Print(const Option_t* delme) const{

    PAngularDistribution::Print();
}


