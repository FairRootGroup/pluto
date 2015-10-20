//TITLE Eta -> pi pi dilepton with form factor and decay angles


 {

     makeDistributionManager();
     
     PReaction my_reaction("_T1=2.2","p","d","He3 eta [e+ e- pi+ pi-]","eta_pi_pi_dilepton");
     my_reaction.Print();
     
     
     //This is for on-line debugging
     TH1F * histo1 = new TH1F ("histo1","dilepton mass",50,0.,0.75);
     
     my_reaction.Do(histo1,"_x = ([e+] + [e-])->M() ");
     
     my_reaction.Loop(300000);
     
     histo1->Draw("");
 }
