//TITLE Demo using the batch syntax to change the omega form factor

 {

     PSimpleVMDFF *ff = new PSimpleVMDFF("vmd_ff_dd@w_to_dilepton_pi0/formfactor",
					 "VMD form factor",-1);
     //use either the AddEquation or the SetVectorMesonMass version
     ff->AddEquation("_ff2 = 0.17918225/( (0.4225- _q2)*(0.4225- _q2) + 0.000676)");
     //ff->SetVectorMesonMass(0.77);


     makeDistributionManager()->Add(ff);


     PReaction my_reaction("2.2","p","p","p p w [dilepton [e+ e-] pi0]");
     my_reaction.Print();

   
     //This is for on-line debugging
     TH1F * histo1 = new TH1F ("histo1","ee invariant mass",100,0.,0.8);
     histo1->Sumw2();
     my_reaction.Do(histo1,"_x = ([e+,1] + [e-,1])->M() ");
     
     my_reaction.Loop(10000);
     
     histo1->Draw("");
     
     
   

 }
