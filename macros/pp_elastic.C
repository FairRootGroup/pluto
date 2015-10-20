//TITLE: pp elastic scattering

{
    gROOT->Reset();


    // pp elastic scattering with beam pz=1.687 GeV/c:

    PParticle *p1=new PParticle("p",0,0,1.687);
    PParticle *p2=new PParticle("p");
    PParticle *p3=new PParticle("p");
    PParticle *p4=new PParticle("p");

    PParticle *q=new PParticle(*p1+*p2);
    PParticle *s[]={q,p3,p4};
    PChannel *c1=new PChannel(s);

    //un-comment next line to have isotropic scattering:
    //makeDistributionManager()->Disable("pp_elastic");

    PChannel *c[]={c1};
    PReaction *r=new PReaction(c,"pp_elastic",1,1);
    r->Print();
    r->loop(500000);

}
