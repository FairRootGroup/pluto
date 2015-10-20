
#include "../../src/PChannelModel.h"

//Template for a formfactor model
//VMD form factor according to Picciotto, Phys.Rev. D45 1569 (1992) / 
//Picciotto and Richardson, Phys.Rev.D48 3395 (1993).


class PEtaPiPiGammaFF : public PChannelModel  {
  
 public:

    using PDistribution::GetWeight;
    PEtaPiPiGammaFF(Char_t *id, Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);
        

 private:

    //Particle parameters:
    Double_t m_pi,m_rho,w_rho;

    ClassDef(PEtaPiPiGammaFF,0)  //Form factor model
};

//---> Do not change something here:
PDistribution* PEtaPiPiGammaFF::Clone(const char*) const {
    //clone the object
    return new PEtaPiPiGammaFF((const PEtaPiPiGammaFF &)* this);
};

PEtaPiPiGammaFF::PEtaPiPiGammaFF(Char_t *id, Char_t *de, Int_t key) : PChannelModel(id, de,key) {
    m_pi = makeStaticData()->GetParticleMass("pi+");
    m_rho = makeStaticData()->GetParticleMass("rho0");
    w_rho = makeStaticData()->GetParticleTotalWidth("rho0");
};

Double_t PEtaPiPiGammaFF::GetWeight(Double_t *mass, Int_t *) {

    Double_t s_pipi = mass[0]*mass[0];

    Double_t Gamma = 0;
    if(s_pipi>=pow(2*m_pi,2))
      Gamma = w_rho*s_pipi/pow(m_rho,2)*pow((1-4*m_pi*m_pi/s_pipi)/(1-4*m_pi*m_pi/m_rho),1.5);

    Double_t FF =  (pow(2*m_rho*m_rho+s_pipi,2)+pow(m_rho*Gamma,2)) / (pow(m_rho*m_rho-s_pipi,2) + pow(m_rho*Gamma,2));

    //cout << s_pipi << endl;
    
    return FF;

}


ClassImp(PEtaPiPiGammaFF)


