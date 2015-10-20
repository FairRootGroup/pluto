//TITLE 

{

//    makeDistributionManager()->Disable("pp_eta_nstar_tcross");
    makeDistributionManager()->Disable("pp_eta_pp_align");
    makeDistributionManager()->Disable("pp_eta_prod_angle");


    PReaction my_reaction(3.953235,"d","p","p p eta (n)");
    //PReaction my_reaction(10.,"d","p","p p eta n");
    my_reaction.Print();   //The "Print()" statement is optional

    //Create my histogram:
    TH1F * histo1 = new TH1F ("histo1","p+p mass",100,2.3,2.8);
    //histo1->Sumw2();
    //Create the container of the histogram list
    PProjector *m1 = new PProjector(); 
    //Dilepton mass
    m1->AddHistogram(histo1,"_x=[p + p]->M()");

    my_reaction.AddBulk(m1);


    my_reaction.Loop(1000);

    histo1->Draw("");

}











