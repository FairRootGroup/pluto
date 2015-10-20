//TITLE  Demonstrate how batch scripting can be used to define filters

{

    //First define our histograms
    TH1F * histo1     = new TH1F ("histo1","p p missing mass",100,0.1,1.0);
    TH1F * histo1_cut = new TH1F ("histo1_cut","p p missing mass",100,0.1,1.0);

    TH1F * histo2      = new TH1F ("histo2","p p missing momentum",100,0.1,1.0);
    TH1F * histo2_cut  = new TH1F ("histo2_cut","p p missing momentum",100,0.1,1.0);
    TH1F * histo2_cut2 = new TH1F ("histo2_cut2","p p missing momentum",100,0.1,1.0);

    TH1F * histo3 = new TH1F ("histo3","Theta of the first proton",100,0,180);

    //Define the reaction
    PReaction my_reaction("3.5","p","p","p p w [pi+ pi- pi0]");
    
    //The batch commands
    my_reaction.Do("pp_miss = [p + p] - ([p,1] + [p,2]); total_miss = [p + p] - ([p,1] + [p,2] + [pi+] + [pi-])");

    //my_reaction.Do("[p,1]->Print();[p,2]->Print(); echo *************");
    //A simple missing mass/momentum histogram
    my_reaction.Do(histo1,"_x= pp_miss->M2()");
    my_reaction.Do(histo2,"_x= pp_miss->P()");

    //Mathematical operations can be used:
    my_reaction.Do(histo3,"_x= ([p,1]->Theta() * 180.)/TMath::Pi()");

    //*******Conditional histogram
    //Only omegas having a certain momentum:
    my_reaction.Do(histo1_cut,"if pp_miss->P() > 0.3; _x= pp_miss->M2()");
    //Control (to be on the save side!):
    my_reaction.Do(histo2_cut,"if pp_miss->P() > 0.3; _x= pp_miss->P()");
    
    
    //*******Conditional Reaction: filters (in our case connected with a control histo)
    //Boolean operations:
    my_reaction.Do(histo2_cut2,"#mom = 0.; if pp_miss->P() > 0.3 && pp_miss->P() < 0.6 ; #mom = 1.");

    my_reaction.Print();

    //Start event loop:
    cout << my_reaction.Loop(50000) << " events recorded" << endl;

    //Plot our on-line histograms
    //******* Starting from here all commands are standard ROOT ***********
    TCanvas *c1 = new TCanvas("c1","Canvas");
    c1->SetLogy(1);
    c1->Divide(2);
    c1->cd(1);
    histo1->Draw("");
    histo1_cut->Draw("same");
    histo1_cut->SetLineColor(2);

    c1->cd(2);
    histo2->Draw("");
    histo2_cut->Draw("same");
    histo2_cut->SetLineColor(2);
    histo2_cut2->Draw("same");
    histo2_cut2->SetLineColor(3);
}
