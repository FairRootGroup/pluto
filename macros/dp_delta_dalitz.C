//TITLE Produce Delta Dalitz events in the quasi-free n+p (d+p) collision
// This macro tests:
//   1. The deuteron Fermi sampling
//   2. The mass-dependent width Breit-Wigner distribution of the Delta 
//   3. The anisotropic (s+p wave) production angle for the pn->pDelta channel
//   4. The dilepton mass for Delta Dalitz decay
{
gROOT->Reset();

PParticle *p1=new PParticle("d",2.);          // projectile   = 1 GeV/u
PParticle *p2=new PParticle("p");             // target
PParticle *p3=new PParticle("p");
PParticle *p4=new PParticle("p");
PParticle *delta=new PParticle("D0");

PParticle *p5=new PParticle("n");
PParticle *dl=new PParticle("dilepton");

PParticle *em=new PParticle("e-");
PParticle *ep=new PParticle("e+");

PParticle *q=new PParticle(*p1+*p2);   // composite p+d
//q->SetW(0.1);

PParticle *s1[]={q,p3,delta,p4}, *s2[]={delta,p5,dl}, *s3[]={dl,em,ep};

PChannel *c1=new PChannel(s1,3);
PChannel *c2=new PChannel(s2);
PChannel *c3=new PChannel(s3);
PChannel *cc[]={c1,c2,c3};

PReaction *r=new PReaction(cc,"dp_delta_dalitz",3,1);

r->Print();
r->loop(10000);
}











