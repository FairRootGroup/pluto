//TITLE A dummy channel model including a new partial decay width

#include "../src/PChannelModel.h"

//Class definition

class PDummyModel : public PChannelModel  {
  
 public:

    PDummyModel(Char_t *id, Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t GetWidth(Double_t mass, Double_t *width, Int_t didx);

 private:

    Double_t parent_pole_mass; //Pole mass of the decay parent
    Int_t parent_id;           //Id of our parent
    Double_t static_width;     //Static width at pole mass of parent
    ClassDef(PDummyModel,0)  //Just a dummy model for mass sampling, NO PHYSICS!
};

PDistribution* PDummyModel::Clone(const char*) const {
    //clone the object
    return new PDummyModel((const PDummyModel &)* this);
};

PDummyModel::PDummyModel(Char_t *id, Char_t *de, Int_t key) : PChannelModel(id, de,key) {
    //Constructor
    if (is_channel<0)
	Warning("PDummyModel","The model (%s) should be bound to CHANNELS only",de);
    
    //Do everything what we need later in advance
    //This saves time in the event loop!

    parent_id = makeStaticData()->GetDecayParent(is_channel);
    parent_pole_mass = makeStaticData()->GetParticleMass(parent_id);
    static_width = makeStaticData()->GetDecayPartialWidth(is_channel);

} ;

Bool_t PDummyModel::GetWidth(Double_t mass, Double_t *width, Int_t didx) {
    //Coupled-channel function for the partial width sampling
    
    //Just a stupid dummy (step) function
    *width = 0;
    if (mass>parent_pole_mass-0.05) *width = static_width;

    return kTRUE;

}


ClassImp(PDummyModel)


