/////////////////////////////////////////////////////////////////////
//Eta double Dalitz decay
//
//Envelope for the complete chain eta->gamma*gamma*->e+e-e+e-
//
//                             Author:  I. Froehlich
//                             Written: 30.3.2010
//                           
//////////////////////////////////////////////////////////////////////

#include "PEtaDoubleDalitzEnv.h"


PDistribution* PEtaDoubleDalitzEnv::Clone(const char*) const {
    //clone the object
    return new PEtaDoubleDalitzEnv((const PEtaDoubleDalitzEnv &)* this);
};

PEtaDoubleDalitzEnv::PEtaDoubleDalitzEnv(const Char_t *id, const Char_t *de, Int_t key) : PChannelModel(id, de,key) {
    //Constructor
    dil1=dil2=parent=NULL;
} ;


Bool_t PEtaDoubleDalitzEnv::Init(void) {
  
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

    ep1  = GetParticle("e+");  
    ep2  = GetParticle("e+");
    em1  = GetParticle("e-");  
    em2  = GetParticle("e-");
  
    return kTRUE;
}


Bool_t PEtaDoubleDalitzEnv::Finalize(void) {
    return kFALSE;
}

Bool_t PEtaDoubleDalitzEnv::EndOfChain(void) {

    //when we are at this stage, all particles are already in lab
    //frame. We habe to boost them into the eta c.m.

    //first make a copy
    PParticle ep1_tmp(ep1);
    PParticle ep2_tmp(ep2);
    PParticle em1_tmp(em1);
    PParticle em2_tmp(em2);

    ep1_tmp.Boost(-parent->BoostVector());  
    ep2_tmp.Boost(-parent->BoostVector());
    em1_tmp.Boost(-parent->BoostVector());  
    em2_tmp.Boost(-parent->BoostVector());

    PParticle dil1_tmp(ep1_tmp);
    dil1_tmp.AddTmp(em1_tmp);
    PParticle dil2_tmp(ep2_tmp);
    dil2_tmp.AddTmp(em2_tmp);

    //rotate such that dileptons point to z-axis
    Double_t Phi = dil1_tmp.Phi();
    Double_t Theta = dil1_tmp.Theta();
    ep1_tmp.RotateZ(-Phi);
    ep1_tmp.RotateY(-Theta);
    dil1_tmp.RotateZ(-Phi);
    dil1_tmp.RotateY(-Theta);
    ep1_tmp.Boost(0,0,-dil1_tmp.Beta());
    //   em1_tmp.RotateZ(-dil1_tmp.Phi());
    //   em1_tmp.RotateY(-dil1_tmp.Theta());
    //   em1_tmp.Boost(0,0,-dil1_tmp.Beta()); --> Unused

    Phi = dil2_tmp.Phi();
    Theta = dil2_tmp.Theta();
    ep2_tmp.RotateZ(-dil2_tmp.Phi());
    ep2_tmp.RotateY(-dil2_tmp.Theta());
    dil2_tmp.RotateZ(-Phi);
    dil2_tmp.RotateY(-Theta);
    ep2_tmp.Boost(0,0,-dil2_tmp.Beta());
    //   em2_tmp.RotateZ(-dil2_tmp.Phi());
    //   em2_tmp.RotateY(-dil2_tmp.Theta());
    //   em2_tmp.Boost(0,0,-dil2_tmp.Beta()); --> Unused

    Double_t theta1=ep1_tmp.Theta();
    Double_t theta2=ep2_tmp.Theta();
    Double_t phi=ep1_tmp.Phi() - ep2_tmp.Phi();

    Double_t sin1_2 = sin(theta1)*sin(theta1);
    Double_t sin2_2 = sin(theta2)*sin(theta2);

    Double_t beta1_2 = 1 - 4* ep1->M2() / dil1_tmp.M2();
    Double_t beta2_2 = 1 - 4* ep2->M2() / dil2_tmp.M2();

    Double_t amplitude =
	0.5*((1 + (1 - beta1_2*sin1_2) *  (1 - beta2_2*sin2_2)) 
	     * sin(phi)* sin(phi) +
	     (2 - beta1_2*sin1_2 - beta2_2*sin2_2) * cos(phi)* cos(phi));

    //cout << "EndOfChain " <<  amplitude  << endl;

    if (PUtils::sampleFlat() > amplitude) return kFALSE;

    return kTRUE;

}

ClassImp(PEtaDoubleDalitzEnv)
