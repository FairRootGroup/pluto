{  // test set of thermal sources   (R.H. 9/8/2000)
    //gROOT->Reset();
    PFireball *source1=new PFireball("pi0",0.,0.1,0,1,0,0.5,0,0,-0.2);
    PFireball *source2=new PFireball("pi+",0.,0.05,0,1,0,0.5,0,0,-0.2);
    PFireball *source3=new PFireball("pi-",0.,0.05,0,1,0,0.5,0,0,-0.2);
    source1->SetW(0.7);  // set weights
    source2->SetW(0.2);
    source3->SetW(0.1);
    source1->Print();
    source2->Print();
    source3->Print();

    PParticle *pi1=new PParticle("pi0");
    PParticle *pi2=new PParticle("pi0");
    PParticle *s1[]={source1,pi1,pi2};
    PParticle *pi3=new PParticle("pi+");
    PParticle *pi4=new PParticle("pi+");
    PParticle *s2[]={source2,pi3,pi4};
    PParticle *pi5=new PParticle("pi-");
    PParticle *pi6=new PParticle("pi-");
    PParticle *pi7=new PParticle("pi-");
    PParticle *s3[]={source3,pi5,pi6,pi7};

    PChannel  *c1=new PChannel(s1,2,1);
    PChannel  *c2=new PChannel(s2,2,1);
    PChannel  *c3=new PChannel(s3,3,1);
    PChannel  *cc[]={c1,c2,c3};
    PReaction *r=new PReaction(cc,"thermal",3,0,0);
    r->Print();
    r->setHGeant(0);   // set to 1, if PLUTO run from HGeant prompt
    r->loop(10000);
}
