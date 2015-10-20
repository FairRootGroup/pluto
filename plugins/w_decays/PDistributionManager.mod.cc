//this line is added to the PDistributionManager.cc


POmega3Pi  *dalitz = new POmega3Pi("w_to_pi+_pi-_pi0_matrix","Omega to 3pi Distribution");

dalitz->Add("w,parent");
dalitz->Add("pi+,daughter");
dalitz->Add("pi-,daughter");
dalitz->Add("pi0,daughter,primary");
dalitz->SetMax(0.2);


makeDistributionManagerUtil()->Add(dalitz);
