// pi- + d -> Delta0 + n -> n + n + pgamma* -> n + n + e- + e+
// This macro tests:
//   1. The beam momentum sampling
//   2. The deuteron Fermi sampling
//   3. The mass-dependent width Breit-Wigner distribution of the Delta 
//   4. The dilepton mass for Delta Dalitz decay
{
gROOT->Reset();

PParticle *p1=new PParticle("pi-",0,0,1.3);   // projectile   = 1.3 GeV/c
PParticle *p2=new PParticle("d");             // target
PParticle *p3=new PParticle("n");
PParticle *delta=new PParticle("D0");

PParticle *p4=new PParticle("n");
PParticle *dl=new PParticle("dilepton");

PParticle *em=new PParticle("e-");
PParticle *ep=new PParticle("e+");

PParticle *q=new PParticle(*p1+*p2);   // composite pi-+d

PParticle *s1[]={q,delta,p3}, *s2[]={delta,p4,dl}, *s3[]={dl,em,ep};

PChannel *c1=new PChannel(s1);
PChannel *c2=new PChannel(s2);
PChannel *c3=new PChannel(s3);
PChannel *cc[]={c1,c2,c3};

PReaction *r=new PReaction(cc,"pid_delta_dalitz",3,1);

r->Print();
r->loop(100000);
}











