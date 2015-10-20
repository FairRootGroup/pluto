{  // test thermal source of pi0's with only Dalitz decays (R.H. 16/6/2005)
    //PData::setFlatMD(1); // use flat Dalitz sampling for better statistics at high mass
    //PFireball *source=new PFireball("pi0",1.0,0.025,0.055,0.5,0.3,0.2,0.,0.,0.);
    PFireball *source=new PFireball("pi0",2.0,0.032,0.089,0.5,0.3,0.2,0.,0.,0.);
    //PFireball *source=new PFireball("pi0",2.0,0.048,0.087,0.9,0.,0.,0.,0.,0.);
    //PFireball *source=new PFireball("pi0",0.0,0.01,0.0,1.0,0.,0.,0.,0.,0.);
    source->Print();
    source->setTrueThermal(kTRUE);
    PParticle *pi0=new PParticle("pi0");
    PParticle *s[]={source,pi0};
    PChannel  *c1=new PChannel(s,1,1,1);
    PParticle *ep = new PParticle("e+");
    PParticle *em = new PParticle("e-");
    PParticle *gam = new PParticle("g");
    PParticle *dilep = new PParticle("dilepton");
    PParticle *pidecay[] = {pi0,gam,dilep};
    PChannel *c2 = new PChannel(pidecay,2,1,1);
    PParticle *dildecay[] = {dilep,ep,em};
    PChannel *c3 = new PChannel(dildecay,2,1,1);
    PChannel  *cc[]={c1,c2,c3};
    PReaction *r=new PReaction(cc,"thermal_pi0",3,0,1);
    r->Print();
    r->setHGeant(0);   // set to 1, if PLUTO run from HGeant prompt
    r->loop(100000);
}
