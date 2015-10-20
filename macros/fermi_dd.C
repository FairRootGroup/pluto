//TITLE <PFermiMomentumDD> The scattering of two nucleons inside deuterons

//This macro demonstrates how the model PFermiMomentumDD
//can be used


#include "../src/PValues.h"

fermi_dd() {

    //This is the real reaction:
    PParticle * beam  =new PParticle("d",2.2);
    PParticle * target=new PParticle("d");

    //Set the values BEFORE using the "+" operator
    //beam->SetValue(P_SCATTER);
    target->SetValue(P_SCATTER);

    PParticle * s= new PParticle(*beam+*target);

    s->Print("scatter");

    //Quasi-free sub-reaction:
    PParticle * beam2  =new PParticle("p");
    PParticle * target2=new PParticle("n");
    PParticle * spectator1=new PParticle("n");
    PParticle * spectator2=new PParticle("p");
    PParticle * s2= new PParticle(*beam2+*target2);

    //The 2 outgoing products of the p-n scattering:
    PParticle * p1  =new PParticle("p");
    PParticle * p2  =new PParticle("n");

    PParticle * cc1[]={s,s2,spectator1,spectator2};
    PParticle * cc2[]={s2,p1,p2};

    PChannel *c1=new PChannel(cc1,3);
    PChannel *c2=new PChannel(cc2,2);
    PChannel *cc[]={c1,c2};

    //PFermiMomentumDD *pmodel = new PFermiMomentumDD("nn_in_dd","Quasi-free particle production");
    PFermiMomentumDD *pmodel = 
	new PFermiMomentumDD("nn_in_dd@d + d_to_p + n_n_p","Quasi-free particle production",-1);

    //now add all particles
    //define spectators and final decay products (the granddaughters)
    pmodel->Add("q,parent");
    pmodel->Add("d,grandparent,beam");
    pmodel->Add("d,grandparent,target");
    pmodel->Add("p,daughter,spectator");
    pmodel->Add("n,daughter,spectator");
    pmodel->Add("q,daughter,composite"); 
    pmodel->Add("p,granddaughter,p1");
    pmodel->Add("n,granddaughter,p2");

    //make it known to the Pluto world:
    makeDistributionManager()->Add(pmodel);
    makeDistributionManager()->Print("user");//The "Print()" statement is optional


    PReaction *r=new PReaction(cc,"fermi_dd",2,1,0,0,0);

    r->Print();


    r->loop(100000);

}
