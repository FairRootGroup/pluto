/////////////////////////////////////////////////////////////////////
//
// Decay Model of Hadron -> 3 unstable hadrons
// 
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PHadronDecayM3.h"


ClassImp(PHadronDecayM3)

PHadronDecayM3::PHadronDecayM3()  {
} ;

PHadronDecayM3::PHadronDecayM3(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de,key) {
    if (is_channel<0)
	Warning("PHadronDecayM3","The model (%s) should be bound to CHANNELS only",de);
  
    //Get particles
    Int_t tid[11];
    tid[0]=10; 
    makeStaticData()->GetDecayModeByKey(primary_key,tid); // retrieve current mode info

    //Parent ALWAYS important (also for the inherited classes)
    parent_id   = makeStaticData()->GetDecayParentByKey(primary_key);
    parent_mass = makeStaticData()->GetParticleMass(parent_id);

    if (tid[0]!=3) 
	Warning("PHadronDecayM3","(%s):  Only 3 body decay",de);


    mass1=makeStaticData()->GetParticleMass(tid[1]);
    mass2=makeStaticData()->GetParticleMass(tid[2]);
    mass3=makeStaticData()->GetParticleMass(tid[3]);
    id1=tid[1];
    id2=tid[2];
    id3=tid[3];
    didx1=-1;
    didx2=-1;
    didx3=-1;
} ;

PDistribution* PHadronDecayM3::Clone(const char*delme) const {
    return new PHadronDecayM3((const PHadronDecayM3 &)* this);
};

Bool_t PHadronDecayM3::Init(void) {
    //Init function called once for each PChannel

    //looking for parent. This is mandatory
    parent = GetParticle("parent");
    if (!parent) {
	Warning("Init","Parent not found");
	return kFALSE;
    }

    //getting 3 daughters
    daughter1 = GetParticle(makeStaticData()->GetParticleName(id1));
    daughter2 = GetParticle(makeStaticData()->GetParticleName(id2));
    daughter3 = GetParticle(makeStaticData()->GetParticleName(id3));


    return kTRUE;
}

Bool_t PHadronDecayM3::Prepare(void) {
    //Things which might change during the eventloop
    
    didx1=daughter1->GetDecayModeIndex(1);
    didx2=daughter2->GetDecayModeIndex(1);
    didx3=daughter3->GetDecayModeIndex(1);

    return kTRUE;


}
 
Double_t PHadronDecayM3::Eval(Double_t x, Double_t y , Double_t z , Double_t t ) const
{
    Double_t res;
    Double_t mass[3];
    mass[0]=x; 
    mass[1]=mass1;
    mass[2]=mass2;
    mass[3]=mass3;

    if (draw_option==0) {
	return ((PChannelModel*)this)->GetWeight(mass);
	//return res;
    }
    if (draw_option==1) {
	((PChannelModel*)this)->GetWidth(x,&res);
	return res;
    }
    if (draw_option==2) {
	((PChannelModel*)this)->GetBR(x,&res);
	return res;
    }
    return 0;
}

int PHadronDecayM3::GetDepth(int i) {
    //check if we have models
    //This also initializes the sub-models

    Int_t a1,a2,a3;
    model1 = makeDynamicData()->GetParticleModel(id1);
    model2 = makeDynamicData()->GetParticleModel(id2);
    model3 = makeDynamicData()->GetParticleModel(id3);
    if (model1) a1=model1->GetDepth(i);
    if (model2) a2=model2->GetDepth(i);
    if (model3) a3=model3->GetDepth(i);

    makeStaticData()->SetDecayEmin(is_channel, mass1+mass2+mass3);

    return TMath::Max(a1+1,TMath::Max(a2+1,a3+1)); 
}

void PHadronDecayM3::SubPrint(Int_t opt) const {
  //Print sub-models
    
    if (model1) {cout << " ";cout << model1->GetDescription();}
    if (model2) {cout << " ";cout << model2->GetDescription();}
    if (model3) {cout << " ";cout << model3->GetDescription();}
}

Bool_t PHadronDecayM3::SampleMass(void) {
    //Mass-sampling wrapper

    Double_t mass[4];
    Int_t didx[4];
    didx[1]=didx1;
    didx[2]=didx2;
    didx[3]=didx3;
    mass[0]=parent->M();
    Bool_t ret = SampleMass(mass,didx);
    daughter1->SetM(mass[1]);
    daughter2->SetM(mass[2]);
    daughter3->SetM(mass[3]);
    return ret;
};

Bool_t PHadronDecayM3::SampleMass(Double_t *mass, Int_t *didx) {
    //Samples the mass of 3 decay products
    //We check if there is a primary model for each
    //decay product, if not, set the mass to the nominal
    //value

    int didx_local1=-1;
    int didx_local2=-1;
    int didx_local3=-1;
    if (didx) didx_local1=didx[1];
    if (didx) didx_local2=didx[2];
    if (didx) didx_local3=didx[3];

    int counter=0;

  repeat1:

    counter++;

    if (model1) model1->SampleMass(&mass[1],&didx_local1);
    else mass[1]=mass1;
    if (model2) model2->SampleMass(&mass[2],&didx_local2);
    else mass[2]=mass2;
    if (model3) model3->SampleMass(&mass[3],&didx_local3);
    else mass[3]=mass3;

    if ((mass[0] < (mass[1] + mass[2] + mass[3])) && (counter <1000)) {
	goto repeat1;
    } else if (counter >=1000){
	return kFALSE;
    }

    //Convolute with 3body phase space
    double m12min2=(mass[1]+mass[2])*(mass[1]+mass[2]);
    double m12max2=(mass[0]-mass[3])*(mass[0]-mass[3]);
    double m23min2=(mass[2]+mass[3])*(mass[2]+mass[3]);
    double m23max2=(mass[0]-mass[1])*(mass[0]-mass[1]);
    double phaseSpace=(m12max2-m12min2)*(m23max2-m23min2);
    m12min2=pow(PData::LMass(id1)+PData::LMass(id2),2);
    m12max2=pow(mass[0]-PData::LMass(id3),2);
    m23min2=pow(PData::LMass(id2)+PData::LMass(id3),2);
    m23max2=pow(mass[0]-PData::LMass(id1),2);
    double phaseSpaceMax=(m12max2-m12min2)*(m23max2-m23min2);
    //    if ( pow(phaseSpace/phaseSpaceMax,0.90)<PUtils::sampleFlat() );
    if ( ((phaseSpace/phaseSpaceMax) <PUtils::sampleFlat() )  && (counter <1000) ) {
	goto repeat1;
    } 	else if (counter >=1000) {
	return kFALSE;
    }
 


    return kTRUE;
};

Bool_t PHadronDecayM3::GetWidth(Double_t mass, Double_t *width, Int_t didx) {
    
    double q_value=mass-(mass1+mass2+mass3);
    double q_value_pole=parent_mass-(mass1+mass2+mass3);
    
    *width=(q_value/q_value_pole)*(q_value/q_value_pole);

    //*width=makeStaticData()->GetDecayPartialWidth(is_channel);
    return kTRUE;

}

Double_t PHadronDecayM3::GetWeight(void) {
    Double_t mass[4];
    mass[0]=parent->M();
    mass[1]=daughter1->M();
    mass[2]=daughter2->M();
    mass[3]=daughter3->M();
    Int_t didx[4];
    didx[1]=didx1;
    didx[2]=didx2;
    didx[3]=didx3;

    //Convolute with 3body phase space
    double m12min2=(mass[1]+mass[2])*(mass[1]+mass[2]);
    double m12max2=(mass[0]-mass[3])*(mass[0]-mass[3]);
    double m23min2=(mass[2]+mass[3])*(mass[2]+mass[3]);
    double m23max2=(mass[0]-mass[1])*(mass[0]-mass[1]);
    double phaseSpace=(m12max2-m12min2)*(m23max2-m23min2);
    m12min2=pow(PData::LMass(id1)+PData::LMass(id2),2);
    m12max2=pow(mass[0]-PData::LMass(id3),2);
    m23min2=pow(PData::LMass(id2)+PData::LMass(id3),2);
    m23max2=pow(mass[0]-PData::LMass(id1),2);
    double phaseSpaceMax=(m12max2-m12min2)*(m23max2-m23min2);

    return GetWeight(mass,didx)*(phaseSpace/phaseSpaceMax); 

}

Double_t PHadronDecayM3::GetWeight(Double_t *mass, Int_t *didx) {
    //Convolution of all 3 masses
    int didx_local1=-1;
    int didx_local2=-1;
    int didx_local3=-1;
    if (didx) didx_local1=didx[1];
    if (didx) didx_local2=didx[2];
    if (didx) didx_local3=didx[3];

    Double_t weight=1.;
    if (model1) weight*=model1->GetWeight(&mass[1],&didx_local1);
    if (model2) weight*=model2->GetWeight(&mass[2],&didx_local2);
    if (model3) weight*=model3->GetWeight(&mass[3],&didx_local3);
    return weight;

}




ClassImp(PHadronDecayM3)
