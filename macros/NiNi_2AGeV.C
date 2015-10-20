{
//
// thermal source model of Ni+Ni 2 AGeV   (R.H. 23/8/2000)
//

Float_t Eb    = 2.0;     // beam energy in AGeV
Float_t T     = 0.089;   // temperature in GeV
Float_t blast = 0.30;    // radial expansion velocity


PFireball *sourceA=new PFireball("p",Eb,T,0.,1.,blast,0.0,0.,0.,0.);
sourceA->Print();
PFireball *sourceB=new PFireball("n",Eb,T,0.,1.,blast,0.0,0.,0.,0.);
sourceB->Print();
PFireball *sourceC=new PFireball("d",Eb,T,0.,1.,blast,0.0,0.,0.,0.);
sourceC->Print();

PFireball *source1=new PFireball("pi0",Eb,T,0.,1.,blast,1.0,0.,0.,0.);
source1->Print();
PFireball *source2=new PFireball("pi+",Eb,T,0.,1.,blast,1.0,0.,0.,0.);
source2->Print();
PFireball *source3=new PFireball("pi-",Eb,T,0.,1.,blast,1.0,0.,0.,0);
source3->Print();
PFireball *source4=new PFireball("eta",Eb,T,0.,1.,blast,1.0,0.,0.,0);
source4->Print();
PFireball *source5=new PFireball("w",Eb,T,0.,1.,blast,1.0,0.,0.,0.);
source5->Print();
PFireball *source6=new PFireball("w",Eb,T,0.,1.,blast,1.0,0.,0.,0.);
source6->Print();
PFireball *source7=new PFireball("rho0",Eb,T,0.,1.,blast,1.0,0.,0.,0.);
source7->Print();
PFireball *source8=new PFireball("phi",Eb,T,0.,1.,blast,1.0,0.,0.,0.);
source8->Print();
PFireball *source9=new PFireball("phi",Eb,T,0.,1.,blast,1.0,0.,0.,0.);
source9->Print();

PFireball *sourceD=new PFireball("D0",Eb,T,0.,1.,blast,0.0,0.,0.,0.);
sourceD->Print();


Float_t Mprot  = 12.3;       // proton
Float_t Mneut  = 13.3;       // neutron
Float_t Mdeut  = 1.2;        // deuteron (= 10% of proton)
                             // <Apart> = 12.3 + 13.3 + 2*1.2 = 28
 
Float_t Mpi0   = 2.77;       // pi0 multiplicity in Ni+Ni 2 AGeV
Float_t Mpip   = 2.63;       // pi+ 
Float_t Mpim   = 2.91;       // pi- 
Float_t Meta   = 0.123;      // eta 
Float_t Momega = 0.015;      // omega 
Float_t Mrho0  = 0.015;      // rho0
Float_t Mphi   = 0.0015;     // phi 
Float_t MDelta = 3./2.*Mpi0; // Delta0 + Delta+
Float_t Mpn    = 1.0e-3;     // pn virtual bremsstrahlung

Float_t enhance = 1.e4;      // enhancement factor

Momega *= enhance;
Mrho0 *= enhance;
Mphi *= enhance;
MDelta *= 100.;
Mpn *= 100.;

sourceA->SetW(1.0);  // set source weights
sourceB->SetW(1.0);
sourceC->SetW(1.0);
source1->SetW(1.0);
source2->SetW(1.0);
source3->SetW(1.0);
source4->SetW(1.0);
source5->SetW(1.0/enhance);
source6->SetW(1.0/enhance);
source7->SetW(1.0/enhance);
source8->SetW(1.0/enhance);
source9->SetW(1.0/enhance);
sourceD->SetW(1.0/100.);
//sourceP->SetW(1.0/enhance);

//                      makeChannel(nMax, nMean, Ebeam)

PChannel *cA = sourceA->makeChannel(30,Mprot);          // p
PChannel *cB = sourceB->makeChannel(30,Mneut);          // n
PChannel *cC = sourceC->makeChannel(5, Mdeut);          // d
PChannel *c1 = source1->makeChannel(9, Mpi0);           // pi0
PChannel *c2 = source2->makeChannel(9, Mpip);           // pi+
PChannel *c3 = source3->makeChannel(9, Mpim);           // pi-
PChannel *c4 = source4->makeChannel(1, Meta);           // eta
PChannel *c5 = source5->makeChannel(1, 5.9e-4*Momega);  // omega-> pi0 e+e-
PChannel *c6 = source6->makeChannel(1, 7.15e-5*Momega); // omega-> e+e-
PChannel *c7 = source7->makeChannel(1, 4.48e-5*Mrho0);  // rho0-> e+e-
PChannel *c8 = source8->makeChannel(1, 3.00e-4*Mphi);   // phi-> e+e-
PChannel *c9 = source9->makeChannel(1, 1.30e-4*Mphi);   // phi-> eta e+e-
PChannel *cD = sourceD->makeChannel(1, 4.4e-5*MDelta);  // Delta-> N e+e-
//PChannel *cP = sourceD->makeChannel(1, Mpn/137.);     // pn-> pn e+e-

PParticle** ptcls;

////////////////    omega Dalitz decay     /////////////////////////////

PParticle* pi01  = new PParticle("pi0");  
PParticle* diel1 = new PParticle("dilepton");
PParticle* elec1 = new PParticle("e-");
PParticle* posi1 = new PParticle("e+");

ptcls = c5->GetParticles();
PParticle* s5_1[] = {ptcls[1], pi01, diel1};
PParticle* s5_2[] = {diel1, elec1, posi1};
PChannel* c5_1 = new PChannel(s5_1,2,1);
PChannel* c5_2 = new PChannel(s5_2,2,1);

////////////////    omega direct decay     /////////////////////////////

PParticle* elec2 = new PParticle("e-");
PParticle* posi2 = new PParticle("e+");

ptcls = c6->GetParticles();
PParticle* s6_1[] = {ptcls[1], elec2, posi2};
PChannel* c6_1 = new PChannel(s6_1,2,1);

////////////////    rho0 direct decay     /////////////////////////////

PParticle* elec3 = new PParticle("e-");
PParticle* posi3 = new PParticle("e+");

ptcls = c7->GetParticles();
PParticle* s7_1[] = {ptcls[1], elec3, posi3};
PChannel* c7_1 = new PChannel(s7_1,2,1);

////////////////    phi direct decay      /////////////////////////////

PParticle* elec4 = new PParticle("e-");
PParticle* posi4 = new PParticle("e+");

ptcls = c8->GetParticles();
PParticle* s8_1[] = {ptcls[1], elec4, posi4};
PChannel* c8_1 = new PChannel(s8_1,2,1);

////////////////    phi Dalitz decay      /////////////////////////////

PParticle* eta1  = new PParticle("eta");
PParticle* diel5 = new PParticle("dilepton");
PParticle* elec5 = new PParticle("e-");
PParticle* posi5 = new PParticle("e+");

ptcls = c9->GetParticles();
PParticle* s9_1[] = {ptcls[1], eta1, diel5};
PParticle* s9_2[] = {diel5, elec5, posi5};
PChannel* c9_1 = new PChannel(s9_1,2,1);
PChannel* c9_2 = new PChannel(s9_2,2,1);

////////////////    Delta Dalitz decay    /////////////////////////////

PParticle* neutron  = new PParticle("n");
PParticle* diel6 = new PParticle("dilepton");
PParticle* elec6 = new PParticle("e-");
PParticle* posi6 = new PParticle("e+");

ptcls = cD->GetParticles();
PParticle* sD_1[] = {ptcls[1], neutron, diel6};
PParticle* sD_2[] = {diel6, elec6, posi6};
PChannel* cD_1 = new PChannel(sD_1,2,1);
PChannel* cD_2 = new PChannel(sD_2,2,1);

//////////////////////////////////////////////////////////////////////

PChannel *cc[] = {
                  cA,cB,cC,           /* nucleons     */
                  c1,c2,c3,c4,        /* pions, eta   */
                  c5,c5_1,c5_2,       /* omega Dalitz */
                  c6,c6_1,            /* omega direct */
                  c7,c7_1,            /* rho0         */
                  c8,c8_1,            /* phi direct   */
                  c9,c9_1,c9_2,       /* phi Dalitz   */
                  cD,cD_1,cD_2        /* Delta Dalitz */
                 };

PReaction *r=new PReaction(cc,"NiNi_2AGeV",22,0,0,0,1);
r->Print();
r->setHGeant(0);   // set to 1, if PLUTO is run from HGeant prompt
r->loop(10000);
}












