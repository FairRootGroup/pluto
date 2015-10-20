/*************************************************************************************
**                                                                                  **
**                                Breakup Reaction:                                 **
**   g+A -> D33+(A-1) -> Delta0+eta+(A-1) -> pi0+N+eta+(A-1) -> 2g+N+2g+(A-1)       **
**             The scattering of a photon and a neutron inside 12C                  **
**      (This macro demonstrates how the model PFermiMomentumGA can be used)        **
**                              *****************                                   **
**                                  Author:                                         **
**                  Lilian Witthauer and Manuel Dieterle                            **
**                                  Date:                                           **
**                                April 2009                                        **
**                             ******************                                   **
**                       (All energies and masses in GeV)                           **
**                                                                                  **
*************************************************************************************/

void C12_breakup_macro(Int_t EventNumber=10000, Double_t StartEBeam = 1.5, Double_t iter=500./1000.) 
{
  char NameFile1[100];
  Double_t EBeamMax = 1.508;                 // Maximum MAMI beam energy
  Double_t EBeam    = StartEBeam;
  Int_t EBeamMeV    = EBeam*1000;
  Double_t DeltaEnergy;

  makeDistributionManager()->Exec("nucleus_fermi:gamma");
 
  Int_t i=0;

  while (EBeam <= EBeamMax)
  {
    // Define projectile (g), target (12C), nucleus fragment (11C)
    PParticle *GammaB  = new PParticle("g",EBeam);
    PParticle *Nucleus = new PParticle("12C");
    PParticle *NucRest = new PParticle("11C");
    
    // IS particle
    PParticle *s = new PParticle(*GammaB+*Nucleus);

    // Quasi-free sub-reaction:
    PParticle *GammaB2  = new PParticle("g");
    PParticle *QFtarget = new PParticle("n"); 
    PParticle *s2       = new PParticle(*GammaB2+*QFtarget);

    // Outgoing products of g-N scattering process:         
    PParticle *Delta   = new PParticle("D0");              
    PParticle *Eta     = new PParticle("eta");             
     
    PParticle *Pi0     = new PParticle("pi0");
    PParticle *Neutron = new PParticle("n");
     
    PParticle *GammaE1   = new PParticle("g");
    PParticle *GammaE2   = new PParticle("g");
    PParticle *GammaP1   = new PParticle("g");
    PParticle *GammaP2   = new PParticle("g");  

    // Define group of particles of each step
    PParticle *cc1[] = {s,NucRest,s2};
    PParticle *cc2[] = {s2,Delta,Eta};
    PParticle *cc3[] = {Delta,Pi0,Neutron};
    PParticle *cc4[] = {Pi0,GammaP1,GammaP2};
    PParticle *cc5[] = {Eta,GammaE1,GammaE2};

    // Allocate groups to decay channel
    PChannel *c1 = new PChannel(cc1,2);
    PChannel *c2 = new PChannel(cc2,2);
    PChannel *ccD0  = new PChannel(cc3,2);
    PChannel *ccPi0 = new PChannel(cc4,2);
    PChannel *ccEta = new PChannel(cc5,2);
    
    // Allocate channels
    PChannel *cc[] = {c1,c2,ccD0,ccPi0,ccEta};
 
    sprintf(NameFile1, "12C_breakup_%i_MeV_in_GeV", EBeamMeV); //File name
    PReaction *Reac = new PReaction(cc,NameFile1,5,1,0,0,0);    // Define reaction
    
    cout << "******************************************" << endl;
    cout << "Beam Energy:" << EBeam << endl;
    cout << "******************************************" << endl; 
    
    Reac->loop(EventNumber);  //Number of events
    Reac->Print();            // Write to .root file 
    

    //Increase energy
    i++;
    DeltaEnergy=i*iter;
    EBeam = StartEBeam+DeltaEnergy;
    EBeamMeV = EBeam*1000;
  }
}
