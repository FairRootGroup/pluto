//TITLE Use the filter file from above

{

    //Please use "make_demo_pfilter.C" before

    makeDistributionManager()->Unpack("pluto_demo_filter.root");

    makeDistributionManager()->Startup("_filter_debug=1");
    makeDistributionManager()->Startup("_filter_exclusive=1");
    //makeDistributionManager()->Startup("_filter_smear_factor=1");

    PReaction my_reaction("2.2","p","p","p p eta [dilepton [e+ e-] g]", "eta_dalitz",1,0,0,0);

    //for the exclusive version:
    TH1F * histo1 = new TH1F ("histo1","pp missing mass",100,0.0,1.);
    my_reaction.Do(histo1,"miss= [p + p]- ( [p,1]+ [p,2] );_x=miss->M()");
    
    
    my_reaction.Print();   //The "Print()" statement is optional
    my_reaction.Loop(10000);

 
}
