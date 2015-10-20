/////////////////////////////////////////////////////////////////////
//Eta double Dalitz decay form factor
//
//References: ???
//
//                             Author:  I. Froehlich
//                             Written: 17.9.2009
//                           
//////////////////////////////////////////////////////////////////////


#include "PEtaDoubleDalitzFF.h"


PDistribution* PEtaDoubleDalitzFF::Clone(const char*) const {
  //clone the object
  return new PEtaDoubleDalitzFF((const PEtaDoubleDalitzFF &)* this);
};

PEtaDoubleDalitzFF::PEtaDoubleDalitzFF(Char_t *id, Char_t *de, Int_t key) : 
  PChannelModel(id, de,key) {
  //Constructor  
  Lambda=0.770; //default value
  dil1=dil2=parent=NULL;
} ;

Bool_t PEtaDoubleDalitzFF::Init(void) {
  
  parent = GetParticle("parent");
  if (!parent) {
    Warning("Init","Parent not found");
    return kFALSE;
  }

  dil1 = GetParticle("dilepton");
  dil2 = GetParticle("dilepton");
  
  if (!dil1 || !dil2) {
    Warning("Init","Dileptons not found");
    return kFALSE;
  }
  
  EnableWeighting();
  SetExpectedWeightMean(-1);

  return kTRUE;
}

Double_t PEtaDoubleDalitzFF::GetWeight() {

  TLorentzVector *P3,*P4;
  P3 = (TLorentzVector*) dil1;
  P4 = (TLorentzVector*) dil2;


  //*******************    Extraction of physics observables      ****************
  Double_t q1_2 =(*P3)*(*P3);//Lorentz Vector of 1st g squared
  Double_t q2_2 =(*P4)*(*P4);//Lorentz Vector of 2nd g squared
  //*************FF from Johan and Fredrik paper***********************
  //     Double_t FF=Lambda^4/(Lambda^2+mVV[0])*Lambda/(Lambda^2+mVV[1]);// including process with Lambda parameter - Ro meson
  //     Double_t FF=Lambda*Lambda*Lambda*Lambda/((Lambda*Lambda-q1_2)*(Lambda*Lambda-q2_2));
  
  //            Double_t FF=Lambda*Lambda/((Lambda*Lambda-q1_2));
  Double_t FF=Lambda*Lambda*Lambda*Lambda/((Lambda*Lambda-q1_2)*(Lambda*Lambda-q2_2));
  float weightFF=FF;

  //weightFF = 20;
  //  cout << weightFF << endl;
  return weightFF;

}


ClassImp(PEtaDoubleDalitzFF)
