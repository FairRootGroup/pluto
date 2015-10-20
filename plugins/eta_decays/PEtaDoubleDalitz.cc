/////////////////////////////////////////////////////////////////////
//Eta double Dalitz decay
//
//References: ???
//
//                             Author:  I. Froehlich
//                             Written: 17.9.2009
//                           
//////////////////////////////////////////////////////////////////////

#include "PEtaDoubleDalitz.h"


PDistribution* PEtaDoubleDalitz::Clone(const char*) const {
  //clone the object
  return new PEtaDoubleDalitz((const PEtaDoubleDalitz &)* this);
};

PEtaDoubleDalitz::PEtaDoubleDalitz(const Char_t *id, const Char_t *de, Int_t key) : PChannelModel(id, de,key) {
  //Constructor
  gRand = new TRandom2();
  dil1=dil2=parent=NULL;
} ;

Double_t PEtaDoubleDalitz::Gen2lepton1(Double_t m)
{
    Double_t MEl=PData::LMass("e+");
    Double_t MEl2=2*MEl;
    Double_t m2=m*m;
    Double_t ARAND[2];
    Double_t WTD;
    Double_t AM;
    do{
        gRand->RndmArray(2,ARAND);
        AM=MEl2*TMath::Exp(TMath::Log(m/MEl2)*ARAND[0]);
        Double_t AM2=AM*AM;
        Double_t RE=MEl2*MEl2/AM2;
        WTD =TMath::Sqrt(1-RE)*(1+RE/2)*TMath::Power(1-AM2/m2,3);
    } while(ARAND[1]>WTD);// ???????????? while(ARAND[1]>WTD)
    return AM;
}


Bool_t PEtaDoubleDalitz::Init(void) {
  
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
    
    formfactor_model=
	GetSecondaryModel("formfactor");
    if (formfactor_model)
	ff_w_max = formfactor_model->GetWeightMax();
    if (ff_w_max<0) {
	Warning("Init","Max value of the FF model not initialized, FF disabled");
	formfactor_model = NULL;	
    } 
    if (!formfactor_model) ff_w_max = 1.;
  
    return kTRUE;
}

Bool_t PEtaDoubleDalitz::SampleMass(void) {

    Double_t MEta = parent->M();
    
    Double_t weight = 1.;
    Double_t mVV[2];
    do {
	mVV[0]=Gen2lepton1(MEta);//IM of the 1-st pair e+e-
	mVV[1]=Gen2lepton1(MEta);//IM of the 2-nd pair e+e-

	if (formfactor_model && formfactor_model->GetVersionFlag(VERSION_SAMPLING)) {
	    weight *= formfactor_model->GetWeight(mVV[0]);
	    weight *= formfactor_model->GetWeight(mVV[1]);
	}

    } while((mVV[0]+mVV[1]>MEta) || ((weight/ff_w_max) < PUtils::sampleFlat())  );
    
    
    dil1->SetM(mVV[0]);
    dil2->SetM(mVV[1]);
    
    return kTRUE;
};

ClassImp(PEtaDoubleDalitz)
