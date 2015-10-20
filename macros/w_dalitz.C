//TITLE Use the long example to produce omega Dalitz events
{
    gROOT->Reset();

    //((PDalitzDecay * )makeDistributionManager()->GetDistribution("w_dalitz"))->SetUseQED(1);
    // use this to get pure QED form factors

    PParticle *pim=new PParticle("pi-",1.1);
    PParticle *p=new PParticle("p");
    PParticle *n=new PParticle("n");
    PParticle *g=new PParticle("dilepton");
    PParticle *w=new PParticle("w");
    PParticle *pi0=new PParticle("pi0");
    PParticle *em=new PParticle("e-");
    PParticle *ep=new PParticle("e+");
    PParticle *q=new PParticle(*pim+*p);
    PParticle *s1[]={q,n,w}, *s2[]={w,pi0,g}, *s3[]={g,em,ep};

    PChannel *ch1=new PChannel(s1);
    PChannel *ch2=new PChannel(s2);
    PChannel *ch3=new PChannel(s3);
    PChannel *ch[]={ch1,ch2,ch3};
    PReaction *r=new PReaction(ch,"w_dalitz",3,1,0,0);



    r->Print();
    // r->SetMaxFileSize(10000);  //make a lot of small files....
    r->loop(10000,1);

    //data.Draw("M()","ID()==51");

}

