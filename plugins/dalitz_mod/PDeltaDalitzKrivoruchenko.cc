/////////////////////////////////////////////////////////////////////
// 
// Dalitz decay following
// Eq. 2 from  Comment on delta radiative and Dalitz decays.
// M.I. Krivoruchenko (Tubingen U. & Moscow, ITEP) , Amand Faessler (Tubingen U.) . Apr 2001. 5pp.
// Published in Phys.Rev.D65:017502,2002.
// e-Print: nucl-th/0104045 (Ref21)
// 
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


#include "PDeltaDalitzKrivoruchenko.h"

PDistribution* PDeltaDalitzKrivoruchenko::Clone(const char*) const {
    //clone the object
    return new PDeltaDalitzKrivoruchenko((const PDeltaDalitzKrivoruchenko &)* this);
};

PDeltaDalitzKrivoruchenko::PDeltaDalitzKrivoruchenko(const Char_t *id, const Char_t *de, Int_t key) : 
    PDalitzDecay(id, de,key) {

    //Constructor
} ;

Bool_t PDeltaDalitzKrivoruchenko:: FreezeOut(void) {
    formfactor_model=
	GetSecondaryModel("formfactor");
    return kTRUE;
}

double PDeltaDalitzKrivoruchenko::dGdM(const int& id, const double& m, const double& ecm) {
    //dGdM from Wolf expression:

    double dgdm = 2.*alpha/(3.*TMath::Pi()*m);

    

    if (m<mass_ee||m>ecm-mass_x) return 0.;
//    cout << m << ":" << ecm << endl;

    double mn=mass_n;
    if (id==Delta_plus || id==S11_plus) 
	mn=mass_p;

    double ff = 2.7*2.7;


    if (formfactor_model)
	ff=formfactor_model->GetWeight(m,ecm);

    dgdm*= (alpha/16.) *
	(ecm + mn) * (ecm + mn) / (ecm*ecm*ecm * mn*mn) *
	sqrt((ecm+mn) * (ecm+mn) - m*m) *
	sqrt(
	    ((ecm-mn) * (ecm-mn) - m*m) *
	    ((ecm-mn) * (ecm-mn) - m*m) *
	    ((ecm-mn) * (ecm-mn) - m*m)
	    ) * ff;

    //   cout << dgdm << endl;

    return dgdm;

}

ClassImp(PDeltaDalitzKrivoruchenko)


