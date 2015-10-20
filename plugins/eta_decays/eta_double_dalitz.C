//TITLE Activate the eta double Dalitz decay eta -> dilepton + dilepton

{

    makeDistributionManager()->Exec("eta_decays");
    makeDistributionManager()->Disable("eta_double_dalitz_complex");
    makeDistributionManager()->Disable("eta_double_dalitz_vmd_ff");

   
    makeDistributionManager()->Print("rare_eta_decays");
    
    PReaction my_reaction("_T1=2.2","p","d","He3 eta[dilepton [e+ e-] dilepton [e+ e-]]","etae+e-e+e-");
    
    
    //This is for on-line tests:
    TH1F * histo1 = new TH1F ("histo1","ee invariant mass",100,0.,0.6);
    TH2F * histo2 = new TH2F ("histo2","ee1 vs. ee2",20,0.,0.6,20,0.,0.6);
    
    my_reaction.Do(histo1,"_x = ([e+,1] + [e-,1])->M() ");
    my_reaction.Do(histo1,"_x = ([e+,2] + [e-,2])->M() ");
    my_reaction.Do(histo2,"_x = ([e+,1] + [e-,1])->M() ; _y =  ([e+,2] + [e-,2])->M()");
    
    
    my_reaction.Print();
    my_reaction.Loop(100000);
    
    
    
    c1 = new TCanvas();
    
    histo1->Draw("");
    
    c2 = new TCanvas();
    
    histo2->Draw("colz");

}



