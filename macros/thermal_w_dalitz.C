//TITLE test thermal source of rho's
// (R.H. 21/3/2006)

{ 

    PFireball *source=new PFireball("w",1.0,0.055,0.0,1.0,0.0,0.,0.,0.,0.); //1AGeV temp.
    source->setTrueThermal(kTRUE);

    PParticle *omegb=new PParticle("w");  // set up omega Dalitz decay
    PParticle *s6[]={source,omegb};
    PChannel  *c14=new PChannel(s6,1,1,1);
    PParticle *ep6 = new PParticle("e+");
    PParticle *em6 = new PParticle("e-");
    PParticle *opi0 = new PParticle("pi0");
    PParticle *dilep4 = new PParticle("dilepton");
    PParticle *omegdecayb[] = {omegb,opi0,dilep4};
    PChannel *c15 = new PChannel(omegdecayb,2,1,1);
    PParticle *dildecay4[] = {dilep4,ep6,em6};
    PChannel *c16 = new PChannel(dildecay4,2,1,1);

    PChannel  *cc[] = {c14,c15,c16};
    PReaction *r = new PReaction(cc,"thermal_w_dalitz",3,1,1);
    r->Print();
    r->SetHGeant(0);   // set to 1, if PLUTO run from HGeant prompt
    r->loop(10000);

    //data.Draw("M()","ID()==52","");

}
