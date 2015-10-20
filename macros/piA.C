//TITLE <PFermiNucleus, PBeamSmearing> Quasi-free pi+A scattering, including beam smearing

{
    
    //Uses PFermiNucleus.C & PCopyBeam
    //This class is here only for didactical reasons. It has no real physics content

    makeStaticData()->AddParticle(90, "A",20.);
    makeStaticData()->AddParticle(91, "A'",19.);
    
    makeStaticData()->AddParticle(90008, "pi+ + A",21.);
    if (!makeStaticData()->IsParticleValid("pi+ + p")) {
	makeStaticData()->AddParticle(14008, "pi+ + p",2.);
    }

    makeStaticData()->AddDecay(-1, "pi+ + A -> (pi+ + p) + A' (quasi-free)", 
			       "pi+ + A","pi+ + p,A'", 1.0 );

    gSystem->CompileMacro( "./PFermiNucleus.C");
    gSystem->CompileMacro( "./PCopyBeam.C");
    
    PFermiNucleus *pmodel = new PFermiNucleus(
			"nucleon_fermi@pi+ + A_to_pi+ + p_A'",
			"Quasi-free particle production",-1);
    
    pmodel->Add("q","parent");
    pmodel->Add("pi+","grandparent","beam");
    pmodel->Add("A","grandparent","target");
    pmodel->Add("A'","daughter","spectator");
    pmodel->Add("q","daughter","composite"); 
    pmodel->Add("p","granddaughter","participant");
    pmodel->Add("pi+","granddaughter","p2");

    makeDistributionManager()->Add(pmodel);
    makeDistributionManager()->Print("user");//The "Print()" statement is optional

    PParticle pi("pi+",2.2); 	//projectil
    PParticle A("A");	        //target 
    PParticle A1("A'");	        //spectator fragment
    
    PParticle piA =pi + A;
    
    PParticle pi1("pi+");
    PParticle p2("p");
    PParticle pip = pi1 + p2;     //quasi-elastic scattering
    
    PParticle p3("p");
    PParticle pi4("pi+");          //outgoing products
   

    PParticle *s0[]={piA,A1,pip};
    PParticle *s1[]={pip,p3,pi4};

    PChannel *c0=new PChannel(s0,2);
    PChannel *c1=new PChannel(s1,2);
    PChannel *cc[]={c0,c1};

    

     PReaction my_reaction(cc,"piA",2,1);

     //Use PCopyBeam macro to copy the original beam + target into the particle array....
     PCopyBeam * copy = new PCopyBeam();
     my_reaction.AddBulk(copy);
     //...could be useful sometimes
     //...this is dangerous if you use GEANT!
  


    makeDistributionManager()->Print("user");//The "Print()" statement is optional
    
    //Create my histograms:
    TH1F * histo2 = new TH1F ("histo2","cos theta",100,-1.,1.);

    //Create the container of the histogram list
    PProjector *m1 = new PProjector(); 
    //cos theta in piA rest frame
    //m1->AddHistogram(histo2,"p1=[p,1]; p1->Boost([pi+ + A]); _x=cos(p1->Theta())");
    m1->AddHistogram(histo2,"p1=[pi+ + p,1]; p1->Boost([pi+ + A]);_x=cos(p1->Theta())");
    m1->Print();

    
    my_reaction.AddBulk(m1);




    my_reaction.Print();


    my_reaction.Loop(1000);


    histo2->Draw();

}
