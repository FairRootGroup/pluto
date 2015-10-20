{  // test thermal source of pi0's   (R.H. 11/8/2000)
//gROOT->Reset();
PFireball *source1=new PFireball("pi0",2.,0.1,0,1,0,0.5,0,0,-0.2);
PFireball *source2=new PFireball("pi+",2.,0.05,0,1,0,0.5,0,0,-0.2);
//source1->SetW(0.90);  // set a weight
source1->setRandomB(40.,40.,0.1);
source2->setRandomB(40.,40.,0.2);
//source1->setMeanN(1.5);
source1->Print();
source2->Print();
PParticle *pi1=new PParticle("pi0");
PParticle *pi2=new PParticle("pi0");
PParticle *pi3=new PParticle("pi0");
PParticle *pi4=new PParticle("pi0");
PParticle *pi5=new PParticle("pi0");
PParticle *pi6=new PParticle("pi0");
PParticle *pi7=new PParticle("pi0");
PParticle *pi8=new PParticle("pi0");
PParticle *pi9=new PParticle("pi0");
PParticle *pi10=new PParticle("pi0");
PParticle *s1[]={source1,pi1,pi2,pi3,pi4,pi5,pi6,pi7,pi8,pi9,pi10};

PParticle *pip1=new PParticle("pi+");
PParticle *pip2=new PParticle("pi+");
PParticle *pip3=new PParticle("pi+");
PParticle *pip4=new PParticle("pi+");
PParticle *pip5=new PParticle("pi+");
PParticle *pip6=new PParticle("pi+");
PParticle *pip7=new PParticle("pi+");
PParticle *pip8=new PParticle("pi+");
PParticle *pip9=new PParticle("pi+");
PParticle *pip10=new PParticle("pi+");
PParticle *pip11=new PParticle("pi+");
PParticle *pip12=new PParticle("pi+");
PParticle *pip13=new PParticle("pi+");
PParticle *pip14=new PParticle("pi+");
PParticle *pip15=new PParticle("pi+");
PParticle *s2[]={source2,pip1,pip2,pip3,pip4,pip5,pip6,pip7,pip8,pip9,pip10,
                 pip11,pip12,pip13,pip14,pip15};
PChannel  *c1=new PChannel(s1,10,1);
PChannel  *c2=new PChannel(s2,15,1);
PChannel  *cc[]={c1,c2};
PReaction *r=new PReaction(cc,"test",2,0,0,0,0);
r->Print();
r->setHGeant(0);   // set to 1, if PLUTO run from HGeant prompt
r->loop(10000,1,0);
}
