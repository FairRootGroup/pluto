{  // test thermal source of Delta's   (R.H. 21/8/2000)
    //gROOT->Reset();
    PFireball *source=new PFireball("D0",1.,0.06,0,1.,0.3,0.5,0,0,0);
    //source->SetW(0.1);  // set a weight
    source->Print();
    source->setMeanN(1.0);
    PParticle *D1=new PParticle("D0");
    PParticle *D2=new PParticle("D0");
    PParticle *D3=new PParticle("D0");
    PParticle *D4=new PParticle("D0");
    PParticle *D5=new PParticle("D0");
    PParticle *D6=new PParticle("D0");
    PParticle *D7=new PParticle("D0");
    PParticle *D8=new PParticle("D0");
    PParticle *D9=new PParticle("D0");
    PParticle *D10=new PParticle("D0");
    PParticle *n1=new PParticle("n");
    PParticle *die1=new PParticle("dilepton");
    PParticle *elec1=new PParticle("e-");
    PParticle *posi1=new PParticle("e+");
    PParticle *n2=new PParticle("n");
    PParticle *die2=new PParticle("dilepton");
    PParticle *elec2=new PParticle("e-");
    PParticle *posi2=new PParticle("e+");
    PParticle *n3=new PParticle("n");
    PParticle *die3=new PParticle("dilepton");
    PParticle *elec3=new PParticle("e-");
    PParticle *posi3=new PParticle("e+");
    PParticle *n4=new PParticle("n");
    PParticle *die4=new PParticle("dilepton");
    PParticle *elec4=new PParticle("e-");
    PParticle *posi4=new PParticle("e+");
    PParticle *n5=new PParticle("n");
    PParticle *die5=new PParticle("dilepton");
    PParticle *elec5=new PParticle("e-");
    PParticle *posi5=new PParticle("e+");

    PParticle *s1[]={source,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10};
    PParticle *s2[]={D1,n1,die1};
    PParticle *s3[]={die1,elec1,posi1};
    PParticle *s4[]={D2,n2,die2};
    PParticle *s5[]={die2,elec2,posi2};
    PParticle *s6[]={D3,n3,die3};
    PParticle *s7[]={die3,elec3,posi3};
    PParticle *s8[]={D4,n4,die4};
    PParticle *s9[]={die4,elec4,posi4};
    PParticle *s10[]={D5,n5,die5};
    PParticle *s11[]={die5,elec5,posi5};
    PChannel *c1=new PChannel(s1,5,1);
    PChannel *c2=new PChannel(s2,2,1);
    PChannel *c3=new PChannel(s3,2,1);
    PChannel *c4=new PChannel(s4,2,1);
    PChannel *c5=new PChannel(s5,2,1);
    PChannel *c6=new PChannel(s6,2,1);
    PChannel *c7=new PChannel(s7,2,1);
    PChannel *c8=new PChannel(s8,2,1);
    PChannel *c9=new PChannel(s9,2,1);
    PChannel *c10=new PChannel(s10,2,1);
    PChannel *c11=new PChannel(s11,2,1);
    //PChannel  *cc[]={c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11};
    PChannel  *cc[]={c1,c2,c4,c6,c8,c10};
    PReaction *r=new PReaction(cc,"thermal_delta",6,0,0,0,0);
    r->Print();
    r->setHGeant(0);   // set to 1, if PLUTO run from HGeant prompt
    r->loop(100000);
}




