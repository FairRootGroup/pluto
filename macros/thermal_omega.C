{  // test thermal source of w's   (R.H. 28/7/2000)

    //omega -> e e decay
    gROOT->Reset();
    PFireball *source=new PFireball("w",2.,0.1,0,1,0,0.5,0,0,-0.2);
    source->Print();
    PParticle *omeg1=new PParticle("w");
    PParticle *s[]={source,omeg1};
    PChannel  *c=new PChannel(s,1,1);
    PChannel  *cc[]={c};
    PReaction *r=new PReaction(cc,"thermal_omega",1,0,0,0,1);
    r->Print();
    r->setHGeant(0);   // set to 1, if PLUTO run from HGeant prompt
    r->loop(50000);
}
