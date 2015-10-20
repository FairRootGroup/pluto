////////////////////////////////////////////////////////
//  
// Adapted TF2
//
////////////////////////////////////////////////////////

#include <math.h>
#include <iostream>
using namespace std;

#ifdef WIN32
#pragma optimize("",off)
#endif

#define INTEGRAL_PRINT_THRESHOLD 10000

#include "PF2.h"
PF2::PF2(const char *name, Double_t (*fcn)(Double_t *, Double_t *), 
	 Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Int_t npar) :
    TF2 (name,fcn,xmin,xmax,ymin,ymax,npar) {
    epsilon=0.000001; //ROOT std.
};


Bool_t PF2::MakeIntegral(void) {
    //Copied from ROOT


    Int_t i,j,cell;
   Double_t dx   = (fXmax-fXmin)/fNpx;
   Double_t dy   = (fYmax-fYmin)/fNpy;
   Int_t ncells = fNpx*fNpy;
 

    if (fIntegral == 0) {
	int print_twentypercent = ncells/5, print_cpc = 1, 
	    print_ipc = print_cpc*print_twentypercent;
	if (ncells>INTEGRAL_PRINT_THRESHOLD) Info("MakeIntegral","Generating array, this can take a while....");

	fIntegral = new Double_t[ncells+1];
	fIntegral[0] = 0;
	Double_t integ;
	Int_t intNegative = 0;
	cell = 0;

	for (j=0;j<fNpy;j++) {
	    for (i=0;i<fNpx;i++) {
		integ = Integral(fXmin+i*dx,fXmin+i*dx+dx,fYmin+j*dy,fYmin+j*dy+dy,epsilon);
		if (cell == print_ipc) {
		    if (ncells>INTEGRAL_PRINT_THRESHOLD) Info("MakeIntegral","...%i%% done",print_cpc*20);
		    print_cpc++;
		    print_ipc = print_cpc*print_twentypercent;
		}
		if (integ < 0) {intNegative++; integ = -integ;}
		fIntegral[cell+1] = fIntegral[cell] + integ;
		cell++;
	    }
	}
	if (intNegative > 0) {
	    Warning("MakeIntegral","function:%s has %d negative values: abs assumed",GetName(),intNegative);
	}
	if (fIntegral[ncells] == 0) {
	    Error("MakeIntegral","Integral of function is zero");
	    return kFALSE;
	}
	for (i=1;i<=ncells;i++) {  // normalize integral to 1
	    fIntegral[i] /= fIntegral[ncells];
	}
	if (ncells>INTEGRAL_PRINT_THRESHOLD) Info("MakeIntegral","...done (%i bins)",ncells);
    }	
    return kTRUE;
};


Double_t PF2::Integral(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t epsilon)
{
    //Copied from ROOT
    // Return Integral of a 2d function in range [ax,bx],[ay,by]
    //
   Double_t a[2], b[2];
   a[0] = ax;
   b[0] = bx;
   a[1] = ay;
   b[1] = by;
   Double_t relerr  = 0;
   Int_t n = 2;
   Int_t minpts = 2*2+2*n*(n+1)+1; //ie 17
   Int_t maxpts = 20*fNpx*fNpy;

   Int_t nfnevl,ifail;

   Double_t result = IntegralMultiple(n,a,b,minpts,maxpts,epsilon,relerr,nfnevl,ifail);

   if (ifail > 0) {
      Warning("Integral","failed code=%d, minpts=%d, maxpts=%d, epsilon=%g, nfnevl=%d, relerr=%g ",ifail,minpts,maxpts,epsilon,nfnevl,relerr);
   }
   return result;
}



ClassImp(PF2)
