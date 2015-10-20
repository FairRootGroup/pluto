/////////////////////////////////////////////////////////////////////
//  Pluto Particle Data
//
//  This class contains tool functions and variables
//
//                             Author:  M.A. Kagarlis
//                             Written: 31/01/99
//                             Revised: 24/05/2004 R.H.
//                             Revised particle tables: 09/06/2004 R.H.
//                             Bug fixes: 13/10/2004 R.H.
//                             mu+mu-: 02/02/2005 R.H.
//                             N1535: 25/01/2007 R.H.
//                             Moved the physics to stand-alone classes 
//                                    (IF, 23.7.2007)
//
//  N.b. that the variables will be moved to the data base soon
//
//  System variables:
//   _system_alloc_verbosity: set =0 to avoid the ALLOCATION printout
//////////////////////////////////////////////////////////////////////

#define INGO_DEBUG

const long double hbar=6.582122e-25; // units of (GeV s)

// local arrays
#include "PData.h"
#include "PKinematics.h"                          // kinematics in-line functions
#include "PUtils.h"                               // utilities class
#include "TF2.h"

#if 0
int PData::useQED=0;        <--- USED???                      // 1 = QED form factors, 0 = VDM (default)
int PData::flatMD=0;        <--- USED???                      // 1 = flat Dalitz mass sampling, 0 = normal
float PData::tempMtScaling=0.;  <--- USED???                  // temperature used for mt scaling
int PData::dbg=0;<--- USED???
int PData::WeightVersion=1;
#endif

int PData::LPW(const int & id, const int & i1, const int & i2) {
    // lowest allowed transition partial wave for the decay
    // of a hadron (id) to two hadrons (i1 & i2).

    int jres=makeStaticData()->GetParticleSpin(id),  // parent 2 x spin
	j1=makeStaticData()->GetParticleSpin(i1),    // 1st particle 2 x spin
	j2=makeStaticData()->GetParticleSpin(i2),    // 2nd particle 2 x spin
	j12max=j1+j2,           // angular momentum limits
	j12min=abs(j1-j2),      // for selection rule
	l=TMath::Min(abs(jres-j12max),abs(jres-j12min))/2,  
	// lowest ang. mom. transfer
	c=(l%2!=abs(makeStaticData()->GetParticleParity(id)
		    -makeStaticData()->GetParticleParity(i1)*
		    makeStaticData()->GetParticleParity(i2))/2);
    // parity correction
    return l+c;
}

int PData::IsMDalitz(const int & idx) {
    // checks if decay mode with index idx is a mesonic Dalitz decay
    // This is calculated here on the fly, dalitz_id is 0 on default

    //TODO: add flag in data base
//    int did=DalitzID[idx];
//    if (did) return (did!=-1) ? did : 0; // -1 means "no"
    
//    DalitzID[idx]=-1;
    if (!makeStaticData()->IsParticleMeson(
	    makeStaticData()->GetDecayParent(idx))) 
	return 0;
    int i, ic[10];
    ic[0]=10; 
    makeStaticData()->GetDecayMode(idx,ic);  // retrieve info for current decay mode
    if (!*ic||*ic>2) return 0; // 0 or >2 products; cannot be Dalitz decay
    for (i=1;i<=2;++i) 
	if (makeStaticData()->IsParticle(ic[i],"dilepton") || 
	    makeStaticData()->IsParticle(ic[i],"dimuon")) {
	//DalitzID[idx]=1;
	return 1;
    }
    return 0;
}
  

#include "TApplication.h"
#include "../Version.h"
#include "../Compiled.h"

#if 0
float PData::impactParam = 0.;
float PData::eventPlane = 0.;
double PData::unstable_width = 0.0001;   <--- USED???
float PData::gVertexX = 0.;
float PData::gVertexY = 0.;
float PData::gVertexZ = 0.;
int   PData::seqNr = -1;
#endif

PSplash *gSplash=0;               // global pointer to PSplash instance
static PSplash PSplashInstance;   // create instance on library load

#include "TLorentzVector.h"
#include "TClass.h"
#include "TMethodCall.h"
#include "TCint.h"

PSplash::PSplash() {
    if (gApplication)
	if (!gApplication->NoLogoOpt()) {
	    cout << "  *********************************************************" << endl;
	    cout << "  * The Pluto event generator                              " << endl;
	    cout << "  * Developed by HADES and all contributing AUTHORS        " << endl;
	    cout << "  * www-hades.gsi.de/computing/pluto/html/PlutoIndex.html  " << endl;
	    cout << "  * Version: " << version_string << endl;
	    cout << "  * Compiled on " << date_string << endl;
	    cout << "  *********************************************************" << endl;

	   	    
	}

}

ClassImp(PSplash)
ClassImp(PData)



    
