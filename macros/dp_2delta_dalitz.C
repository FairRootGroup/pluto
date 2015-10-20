// d + p -> Delta0 + Delta+ -> p + n + p + 2pgamma* -> p + n + p + 2e- + 2e+
// This macro tests:
//   1. The deuteron Fermi sampling
//   2. The mass-dependent width Breit-Wigner distribution of the Delta 
//   3. The anisotropic (s+p wave) production angle for the pn->pDelta channel
//   4. The dilepton mass for Delta Dalitz decay
{
gROOT->Reset();
PParticle *p1=new PParticle("d",4.);          // projectile   = 2 GeV/u
PParticle *p2=new PParticle("p");             // target
PParticle *p3=new PParticle("p");
PParticle *p4=new PParticle("p"); 

PParticle *delta=new PParticle("D0");
PParticle *delta2=new PParticle("D+");

PParticle *p5=new PParticle("n");
PParticle *dl=new PParticle("dilepton");
PParticle *dl2=new PParticle("dilepton");

PParticle *em=new PParticle("e-");
PParticle *ep=new PParticle("e+"); 
PParticle *em2=new PParticle("e-");
PParticle *ep2=new PParticle("e+"); 

PParticle *q=new PParticle(*p1+*p2);   // composite p+d

PParticle *s1[]={q,delta,delta2,p4}, *s2[]={delta,p5,dl}, *s3[]={dl,em,ep},
          *s4[]={delta2,p3,dl2}, *s5[]={dl2,em2,ep2};

PChannel *c1=new PChannel(s1,3);
PChannel *c2=new PChannel(s2);
PChannel *c3=new PChannel(s3);
PChannel *c4=new PChannel(s4);
PChannel *c5=new PChannel(s5);
PChannel *cc[]={c1,c2,c3,c4,c5};

PReaction *r=new PReaction(cc,"dp_2delta_dalitz",5,1);

r->Print();
r->loop(100);
}











