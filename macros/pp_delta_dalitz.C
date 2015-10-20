// p + p -> p + Delta+ -> p + p + gamma* -> p + p + e- + e+
// This macro tests:
//   1. The mass-dependent width Breit-Wigner distribution of the Delta 
//   2. The anisotropic (s+p wave) production angle for the pp->pDelta channel
//   3. The dilepton mass for Delta Dalitz decay
{
    gROOT->Reset();
    
    // use to get pure QED form factors:
    ((PDalitzDecay * )makeDistributionManager()->GetDistribution("D+_dalitz"))->SetUseQED(1);
    
    
    
    PParticle *p1=new PParticle("p",2.2);
    PParticle *p2=new PParticle("p"),*p3=new PParticle("p"),*p4=new PParticle("p");
    PParticle *delta=new PParticle("D+");
    PParticle *g=new PParticle("dilepton");
    PParticle *em=new PParticle("e-");
    PParticle *ep=new PParticle("e+");
    PParticle *q=new PParticle(*p1+*p2);
    PParticle *s1[]={q,p3,delta}, *s2[]={delta,p4,g}, *s3[]={g,em,ep};
    PChannel *c1=new PChannel(s1), *c2=new PChannel(s2), *c3=new PChannel(s3), *cc[]={c1,c2,c3};
    PReaction *r=new PReaction(cc,"pp_delta_dalitz",3,1);
    
    makeDistributionManager()->GetDistribution("D+_dalitz")->Print();
    
    r->Print();
    r->loop(100000);
}
