//TITLE Eta -> pi pi gamma with matrix element


 {

     //makeDistributionManager()->Disable("eta_decays");
     //makeDistributionManager()->Disable("eta_pipi_gamma_matrix_weighting");
     //One can enable weighting as an option, this can be faster...
     //makeDistributionManager()->Enable("eta_pipi_gamma_matrix_weighting");
     //...but one has to use the W() of the PParticle later
     
     makeDistributionManager();

     //Compile and attach FF model:
     gSystem->CompileMacro( "./PEtaPiPiGammaFF.C");
     PEtaPiPiGammaFF * newmodel = new PEtaPiPiGammaFF("my_ff@eta_to_pi+_pi-_gamma/formfactor",
						      "My own FF for eta -> pi+ pi- gamma",-1);
     newmodel->SetWeightMax(10.);
     makeDistributionManager()->Add(newmodel);

     //for debugging:
     //newmodel->SetRange(0.1,0.75);
     //newmodel->Draw();


     
     PReaction my_reaction(2.2,"p","p","p p eta [g pi+ pi-]","eta_pi_pi_gamma");
     my_reaction.Print();
     
     
     //This is for on-line debugging
     TH2F * histo2 = new TH2F ("histo2","Dalitz",20,0.,0.2,20,0.,0.2);
     TH1F * histo1 = new TH1F ("histo1","pipi mass",50,0.25,0.75);
     
     my_reaction.Do(histo2,"_x = ([pi+,1] + [g,1])->M2()  ; _y = ([pi-,1] + [g,1])->M2()");
     my_reaction.Do(histo1,"_x = ([pi+,1] + [pi-,1])->M() ");
     
     my_reaction.Loop(10000);
     
     
     //     histo2->Draw("box");
     histo1->Draw("");
 }
