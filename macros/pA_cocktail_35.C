//TITLE Thermal source model of p+nucleus 3.5 AGeV 

{
// To save space, only leptons branches are generated and written to file.
// allows for different seeds each time the macro is started

PUtils::SetSeed(0);

// makeDistributionManager()->Disable("helicity_angles");

Bool_t freeze=kFALSE;

Float_t Eb    = 3.5;     // beam energy in AGeV
Float_t ctr   = 1.0;     // 
//Float_t bmax  = 3.9;     // max. impact parameter (corresponds to 60% of all)
Float_t bmax  =0.0;         // Poisson sampling
if (bmax>0.) ctr = 1.;   // use b sampling instead of Poisson sampling

Float_t T1    = 0.080;   // temperature in GeV (for pion 2-component spectra)
Float_t T2    = 0.080;   // temperature in GeV (assume this for thermalized source)
Float_t frac  = 1.0;     // fraction of pion low-T component (from Jehad's fit to QMD)
Float_t blast = 0.0;     // radial expansion velocity

Float_t A2    = 1.0;     // polar distribution (from KaoS pion data)
Float_t A4    = 0.0;
Float_t v1    = 0.0;     // side flow
Float_t v2    = 0.0;     // elliptic flow

// Multiplicties for 3.5 AGeV P+nucleus min. bias events
Float_t Mprot  = 2.03*ctr;       // proton
Float_t Mneut  = 2.03*ctr;       // neutron
Float_t Mdeut  = 0.0*ctr;      // deuteron (= 10% of proton)
                                // <Apart> = ctr*(8.9 + 8.9 + 2*0.89) = 19.6 = A/2
Float_t Mpi0   = 0.60*ctr;      // pi0 multiplicity in p+C 3.5 AGeV
Float_t Mpip   = 0.60*ctr;      // pi+
Float_t Mpim   = 0.60*ctr;      // pi-
Float_t Meta   = 0.031*ctr;     // eta
Float_t Momega = 0.011*ctr;     // omega
Float_t Mrho0  = 0.011*ctr;     // rho0  
Float_t Mphi   = 0.0005*ctr;    // phi
Float_t MDelta = 2.*3./2.*Mpi0; // Delta0 + Delta+  

Float_t enhance = 500.;           // enhancement factor for mesons
//Float_t enhance = 5000.; 
// enhance = 1.;           // enhancement factor
Float_t enhancepi=4.;            //enhancement factor for pi
//Float_t enhancepi=1;            //enhancement factor for pi

Mpi0 *=enhancepi;
//Mpip *=enhancepi;
Meta *=enhance;
Momega *= enhance;
Mrho0 *= enhance;
Mphi *= enhance;
MDelta *= enhance;

PFireball *source1=new PFireball("pi0",Eb,T1,T2,frac,blast,A2,A4,v1,v2);
PFireball *source2=new PFireball("eta",Eb,T2,0.,1.,blast,0.0,A4,v1,v2);
PFireball *source3=new PFireball("w",Eb,T2,0.,1.,blast,0.0,A4,v1,v2);
PFireball *source4=new PFireball("w",Eb,T2,0.,1.,blast,0.0,A4,v1,v2);

// source4->getFuncE()->Draw();
// return;

PFireball *source5=new PFireball("rho0",Eb,T2,0.,1.,blast,0.0,A4,v1,v2);
PFireball *source6=new PFireball("phi",Eb,T2,0.,1.,blast,0.0,A4,v1,v2);
PFireball *source7=new PFireball("phi",Eb,T2,0.,1.,blast,0.0,A4,v1,v2);
PFireball *source8=new PFireball("D0",Eb,T2,0.,1.,blast,0.0,A4,v1,v2);
PFireball *source9=new PFireball("phi",Eb,T1,T2,frac,blast,A2,A4,v1,v2);
PFireball *source10=new PFireball("pi-",Eb,T1,T2,frac,blast,A2,A4,v1,v2);
PFireball *source11=new PFireball("p",Eb,T1,T2,frac,blast,A2,A4,v1,v2);
PFireball *source12=new PFireball("eta",Eb,T2,0.,1.,blast,0.0,A4,v1,v2);

source1->SetW(1.0/enhancepi);
source2->SetW(1.0/enhance);
source3->SetW(1.0/enhance);
source4->SetW(1.0/enhance);
source5->SetW(1.0/enhance);
source6->SetW(1.0/enhance);
source7->SetW(1.0/enhance);
source8->SetW(1.0/enhance);
source9->SetW(1.0);
source10->SetW(1.0);
source11->SetW(1.0);
source12->SetW(1.0/enhance);


PChannel *c1 = source1->makeChannel(20, 0.012*Mpi0);          // pi0 decay directly

PChannel *c2 = source2->makeChannel(1, 0.006*Meta);     // eta -> gamma e+e-
PChannel *c3 = source3->makeChannel(1, 5.9e-4*Momega);  // omega -> pi0 e+e-
PChannel *c4 = source4->makeChannel(1, 7.15e-5*Momega); // omega -> e+e-
PChannel *c5 = source5->makeChannel(1, 4.48e-5*Mrho0);  // rho0 -> e+e-
PChannel *c6 = source6->makeChannel(1, 3.00e-4*Mphi);   // phi -> e+e-
PChannel *c7 = source7->makeChannel(1, 1.30e-4*Mphi);   // phi -> eta e+e-
PChannel *c8 = source8->makeChannel(1, 4.4e-5*MDelta);  // Delta -> N e+e-
PChannel *c9 = source9->makeChannel(1,0.49*Mphi); //phi->K+K-
/*
PChannel *c9 = source9->makeChannel(10, Mpip);  // Pion plus
PChannel *c10 = source10->makeChannel(10, Mpim);  // Pion minus
PChannel *c11 = source11->makeChannel(10, Mprot);  // Proton
*/
PChannel *c12 = source12->makeChannel(1, 0.21*Meta);     // eta -> gamma, gamma 
//PChannel *c14 = source13->makeChannel(1, 0.99*Mpi0);     // pi0 -> gamma, gamma
//

source1->Print();
source2->Print();
source3->Print();
source4->Print();
source5->Print();
source6->Print();
source7->Print();
source8->Print();
source9->Print();
source10->Print();
source11->Print();
source12->Print();


PParticle** ptcls;

/// Dalitz decays

PParticle* gam1  = new PParticle("g");
PParticle* diel1 = new PParticle("dilepton");
PParticle* elec1 = new PParticle("e-");
PParticle* posi1 = new PParticle("e+");

////////////////    pi0 Dalitz decay       /////////////////////////////

ptcls = c1->GetParticles();
PParticle* s1_1[] = {ptcls[1], gam1, diel1};
PParticle* s1_2[] = {diel1, elec1, posi1};
PChannel* c1_1 = new PChannel(s1_1,2,1);
PChannel* c1_2 = new PChannel(s1_2,2,1);

/// Dalitz decays

PParticle* gam2  = new PParticle("g");
PParticle* diel2 = new PParticle("dilepton");
PParticle* elec2 = new PParticle("e-");
PParticle* posi2 = new PParticle("e+");

////////////////    eta Dalitz decay       /////////////////////////////

ptcls = c2->GetParticles();
PParticle* s2_1[] = {ptcls[1], gam2, diel2};
PParticle* s2_2[] = {diel2, elec2, posi2};
PChannel* c2_1 = new PChannel(s2_1,2,1);
PChannel* c2_2 = new PChannel(s2_2,2,1);


////////////////    omega Dalitz decay     /////////////////////////////

PParticle* pi03  = new PParticle("pi0");
PParticle* diel3 = new PParticle("dilepton");
PParticle* elec3 = new PParticle("e-");
PParticle* posi3 = new PParticle("e+");

ptcls = c3->GetParticles();
PParticle* s3_1[] = {ptcls[1], pi03, diel3};
PParticle* s3_2[] = {diel3, elec3, posi3};
PChannel* c3_1 = new PChannel(s3_1,2,1);
PChannel* c3_2 = new PChannel(s3_2,2,1);

////////////////    omega direct decay     /////////////////////////////

PParticle* elec4 = new PParticle("e-");
PParticle* posi4 = new PParticle("e+");

ptcls = c4->GetParticles();
PParticle* s4_1[] = {ptcls[1], elec4, posi4};
PChannel* c4_1 = new PChannel(s4_1,2,1);

////////////////    rho0 direct decay     /////////////////////////////

PParticle* elec5 = new PParticle("e-");
PParticle* posi5 = new PParticle("e+");

ptcls = c5->GetParticles();
PParticle* s5_1[] = {ptcls[1], elec5, posi5};
PChannel* c5_1 = new PChannel(s5_1,2,1);

////////////////    phi direct decay      /////////////////////////////

PParticle* elec6 = new PParticle("e-");
PParticle* posi6 = new PParticle("e+");

ptcls = c6->GetParticles();
PParticle* s6_1[] = {ptcls[1], elec6, posi6};
PChannel* c6_1 = new PChannel(s6_1,2,1);

////////////////    phi Dalitz decay      /////////////////////////////

PParticle* eta7  = new PParticle("eta");
PParticle* diel7 = new PParticle("dilepton");
PParticle* elec7 = new PParticle("e-");
PParticle* posi7 = new PParticle("e+");

ptcls = c7->GetParticles();
PParticle* s7_1[] = {ptcls[1], eta7, diel7};
PParticle* s7_2[] = {diel7, elec7, posi7};
PChannel* c7_1 = new PChannel(s7_1,2,1);
PChannel* c7_2 = new PChannel(s7_2,2,1);

////////////////    Delta Dalitz decay    /////////////////////////////

PParticle* neut8 = new PParticle("n");
PParticle* diel8 = new PParticle("dilepton");
PParticle* elec8 = new PParticle("e-");
PParticle* posi8 = new PParticle("e+");

ptcls = c8->GetParticles();
PParticle* s8_1[] = {ptcls[1], neut8, diel8};
PParticle* s8_2[] = {diel8, elec8, posi8};
PChannel* c8_1 = new PChannel(s8_1,2,1);
PChannel* c8_2 = new PChannel(s8_2,2,1);

//phi->K+,k- decay

ptcls= c9->GetParticles();
PParticle* kp    = new PParticle("K+");
PParticle* km    = new PParticle("K-");
PParticle* s9_1[] = {ptcls[1],kp,km};
PChannel* c9_1 = new PChannel(s9_1,2,1);

//// eta two photon decay   ///////////////////////////////////////

PParticle* gamma1 = new PParticle("g");
PParticle* gamma2 = new PParticle("g");

ptcls = c12->GetParticles();
PParticle* s12_1[] = {ptcls[1], gamma1, gamma2};
PChannel* c12_1 = new PChannel(s12_1,2,1);

//////////////////////////////////////////////////////////////////////




PChannel *cc[] = {
                  c1,c1_1,c1_2,        /* pi0          */
                  c2,c2_1,c2_2,       /* eta Dalitz   */
                  c3,c3_1,c3_2,       /* omega Dalitz */
                  c4,c4_1,            /* omega direct */
                  c5,c5_1,            /* rho0         */
                  c6,c6_1,            /* phi direct   */
                  c7,c7_1,c7_2,       /* phi Dalitz   */
                  c8,c8_1,c8_2,        /* Delta Dalitz */
//                 c9,c9_1,                  /* phi K+,k- production */
//                  c10,                 /*Pim production*/
//                  c11,                 /* Proton production */
//                  c12,c12_1            /* eta 2 photon decay */
                 };
//change  *cc[]if you add want to add new channels and modify param after outFile
 
PReaction *r=new PReaction(cc,"pA",21,1,0,0,1);
//PReaction *r=new PReaction(cc,"pA",6,1,0,0,1);
 r->Print();

//r->setUserSelection(selectLeptons);          // this works for a loaded function
//r->setUserSelection(&select);                  //calls filter function in PReaction.cc
//r->setTrigCond(1);    // require at least 1 lepton per event

//r->SetDecayAll(1.);   // decay all particles with tau < 1 ns
r->SetHGeant(0);      // set to 1, if PLUTO is run from HGeant prompt
r->loop(10000);   // do n events

//data.Draw("M()","(ID() == 51 || ID()==52 || ID()==41) * W()");
//data.Draw("M()","(ID()==41) * W()","same"); 

}

