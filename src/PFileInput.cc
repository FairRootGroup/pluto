////////////////////////////////////////////////////////
//  File input interface implementation file
//
//  This class serves as interface to event input from ascii file, e.g.
//  HGeant, UrQMD, etc. files.
//  The constructor opens the file in read-only mode.
//  The PChannel::decay() function reads then one event and fills ptcls[]
//
//                    Author:  Romain Holzmann
//                    Written: 08/11/2002
//                    Revised: 08/11/2002 RH
//                    Revised: 22/03/2005 RH
//
////////////////////////////////////////////////////////

#include "PUtils.h"
#include "PFileInput.h"

ClassImp(PFileInput)

PFileInput:: PFileInput(Char_t* mode, Char_t* filename) :
                                                         PParticle("File") {
// 
// presently mode can be "HGeant" 
//
  if (strcmp(mode,"HGeant")==0) imode = 1;      // HGeant event input format 
  else if (strcmp(mode,"UrQMD")==0) imode = 2;  // UrQMD format (not impl.)
  else if (strcmp(mode,"HSM")==0) imode = 3;    // HSM format (not impl.)
  else imode = 0;
  name=mode;
  part=NULL;
  fp = fopen(filename,"r");   // open event file

  if (imode==0 || fp==NULL)
                  printf("\n*** Error in PFileInput constructor ***\n");
  SetPxPyPzE(0.,0.,0.,1.);
  wt = 1.;
  file = filename;
}

void PFileInput::setToMidrapidity(float agev) {  // cm->lab transformation
  // midrapidity with agev the energy/nucleon of the projectile (in GeV/u)
   if (agev <= 0.0) return;
   double bx = 0.;
   double by = 0.;
   double bz = sqrt(agev/(agev+2.*0.9315));
   Boost(bx,by,bz);
}

Int_t PFileInput::readEventHeader(Float_t& b) {  // read event header
   Int_t ret, evNb, nPar;
   Float_t Ebeam;
   switch (imode) {
   case 1:
     ret = fscanf(fp,"%i %i %f %f %i", &evNb, &nPar, &Ebeam, &b, &flag);
     if (ret==EOF) return -1;
     else return nPar;
     break;
   case 2:  // not implemented
     break;
   case 3:  // not implemented
     break;
   default: ;
   }
   return -1;
}

Int_t PFileInput::readParticle(Double_t& px, Double_t& py, Double_t& pz,
       Double_t &E, Double_t& vx, Double_t& vy, Double_t& vz, Double_t& vt,
       Int_t& Id, Int_t& srcId, Int_t& parId, Int_t& parInd, Double_t& weight) {
   Int_t ret = -1;
   switch (imode) {
   case 1:
     vx = vy = vz = vt = 0.;
     // flag decides how much info is read per particle
     if (flag==0) ret = fscanf(fp,"%le %le %le %le %i", &E, &px, &py, &pz, &Id);
     else if (flag==1) ret = fscanf(fp,"%le %le %le %le %i %le",
		       &E, &px, &py, &pz, &Id, &weight); // + weight
     else if (flag==2) ret = fscanf(fp,"%le %le %le %le %i %i %i %le",
		       &E, &px, &py, &pz, &Id, &srcId, &parId, &weight); // + ids
     else if (flag==-2) ret = fscanf(fp,"%le %le %le %le %i %i %i %i %le",
		        &E, &px, &py, &pz, &Id, &srcId, &parId, &parInd, &weight); // + ids
     else if (flag==3) ret = fscanf(fp,"%le %le %le %le %le %le %le %i %i %i %le",
		       &E, &px, &py, &pz, &vx, &vy, &vz, &Id, &srcId, &parId, &weight); // + vertex
     else if (flag==-3) ret = fscanf(fp,"%le %le %le %le %le %le %le %i %i %i %i %le",
		        &E, &px, &py, &pz, &vx, &vy, &vz, &Id, &srcId, &parId, &parInd, &weight); // + vertex
     else if (flag==4) ret = fscanf(fp,"%le %le %le %le %le %le %le %le %i %i %i %le",
		       &E, &px, &py, &pz, &vx, &vy, &vz, &vt, &Id, &srcId, &parId, &weight); // + time
     else if (flag==-4) ret = fscanf(fp,"%le %le %le %le %le %le %le %le %i %i %i %i %le",
		        &E, &px, &py, &pz, &vx, &vy, &vz, &vt, &Id, &srcId, &parId, &parInd, &weight); // + time
     break;
   case 2:  // not implemented
     break;
   case 3:  // not implemented
     break;
   default: ;
   }
   if (E<0.0) E = -E;
   if (ret==EOF) ret = -1;
   return ret;
}

void PFileInput:: Print(const Option_t* delme) const {
  printf(" Input type: %s\n",name);
  printf(" File name: %s\n",file);
  printf(" rapidity=%f  weight=%f\n", Rapidity(), W());
}

PChannel* PFileInput::makeChannel(Int_t nMax, Float_t Ebeam) {
  //
  // set up a reaction channel for this interface
  //
   if (Ebeam>0.) setToMidrapidity(Ebeam);

   part = new PParticle*[nMax+1];
   part[0] = this;
   for (Int_t i=1;i<=nMax;i++) {
     part[i] = new PParticle("dummy"); 
   }
   PChannel* chan = new PChannel(part,nMax,1);
   return chan;
}




