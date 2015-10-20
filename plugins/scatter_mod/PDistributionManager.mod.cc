//this line is added to the PDistributionManager.cc

{
    //offer new pp distribution as an alternative
    TF1 * dummy = new TF1("dummy","1",-1,1);
    gROOT->GetListOfFunctions()->Remove(dummy);
    PSaidLowEnergy * newsaid = new PSaidLowEnergy("low_energy_pp_elastic",
						  "Low energy scattering <scatter_mod>");
    newsaid->Add("q"   , "PARENT");
    newsaid->Add("p"   , "GRANDPARENT", "beam");
    newsaid->Add("p"   , "GRANDPARENT", "target");
    newsaid->Add("p, daughter");
    newsaid->Add("p, daughter, primary");
    newsaid->SetAngleFunction(dummy);
    makeDistributionManagerUtil()->Add(newsaid);
    makeDistributionManagerUtil()->AlternativeTo("low_energy_pp_elastic","pp_elastic");

    PSaidPN * pnsaid = new PSaidPN("pn_elastic",
				   "p+n scattering <scatter_mod>");
    pnsaid->Add("q"   , "PARENT");
    pnsaid->Add("p"   , "GRANDPARENT", "beam");
    pnsaid->Add("n"   , "GRANDPARENT", "target");
    pnsaid->Add("n, daughter");
    pnsaid->Add("p, daughter, primary");
    pnsaid->SetAngleFunction(dummy);
    makeDistributionManagerUtil()->Add(pnsaid);


}

//PluginInfo("...");

