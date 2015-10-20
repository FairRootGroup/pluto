
{

    //Uses PFermiNucleus.C
    //This class is here only for didactical reasons. It has no real physics content
    //See nucleus_fermi plugin for a real implementation

    makeStaticData()->AddParticle(90, "A",20.);
    makeStaticData()->AddParticle(91, "A'",19.);
    
    makeStaticData()->AddParticle(90014, "p + A",21.);
    if (!makeStaticData()->IsParticleValid("p + p")) {
	makeStaticData()->AddParticle(14014, "p + p",2.);
    }

    makeStaticData()->AddDecay(-1, "p + A -> (p + p) + A' (quasi-free)", 
			       "p + A","p + p,A'", 1.0 );

    gSystem->CompileMacro( "./PFermiNucleus.C");
    
    PFermiNucleus *pmodel = new PFermiNucleus(
			"n_fermi@p + A_to_p + p_A'",
			"Quasi-free particle production",-1);
    
    pmodel->Add("q","parent");
    pmodel->Add("p","grandparent","beam");
    pmodel->Add("A","grandparent","target");
    pmodel->Add("A'","daughter","spectator");
    pmodel->Add("q","daughter","composite"); 
    pmodel->Add("p","granddaughter","participant");
    pmodel->Add("p","granddaughter","p2");

    makeDistributionManager()->Add(pmodel);
    makeDistributionManager()->Print("user");//The "Print()" statement is optional

    makeDistributionManager()->Disable("pp_elastic");

    PParticle p("p",2.2); 	//projectil
    PParticle A("A");	        //target 
    PParticle A1("A'");	        //spectator fragment
    
    PParticle pA =p + A;
    
    PParticle p1("p");
    PParticle p2("p");
    PParticle pp = p1 + p2;     //quasi-elastic scattering
    
    PParticle p3("p");
    PParticle p4("p");          //outgoing products
   

    PParticle *s0[]={pA,A1,pp};
    PParticle *s1[]={pp,p3,p4};

    PChannel *c0=new PChannel(s0,2);
    PChannel *c1=new PChannel(s1,2);
    PChannel *cc[]={c0,c1};


    PReaction my_reaction(cc,"pA",2,1);

  
    makeDistributionManager()->Print("user");//The "Print()" statement is optional
    
    //Create my histograms:
    TH1F * histo1 = new TH1F ("histo1","pp mass",100,2.5,4.0);
    TH1F * histo2 = new TH1F ("histo2","cos theta",100,-1.,1.);
    histo1->Sumw2();
    //Create the container of the histogram list
    PProjector *m1 = new PProjector(); 
    //Dilepton mass
    m1->AddHistogram(histo1,"_x=[p + p]->M()");
    //cos theta in pA rest frame
    m1->AddHistogram(histo2,"p1=[p,1]; p1->Boost([p + A]); _x=cos(p1->Theta())");


    
    my_reaction.AddBulk(m1);




    my_reaction.Print();


    my_reaction.Loop(10000);

    TCanvas *can1 = new TCanvas();
    
    histo1->Draw();
    
    TCanvas *can2 = new TCanvas();
    histo2->Draw();

}
