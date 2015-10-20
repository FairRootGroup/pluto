//TITLE Low-energy pp scattering

{
  gROOT->Reset();

  makeDistributionManager()->Enable("low_energy_pp_elastic");
  
  TCanvas *test_elastic_phi = new
    TCanvas ("test_elastic_phi","lowE pp elastic",100,10,600,400);
  
  TCanvas *test_elastic_theta = new
    TCanvas ("test_elastic_theta","lowE pp elastic",200,5,700,800);
  test_elastic_theta->Divide(2,2);
  
  PReaction ela_pp("0.05","p","p","p p");
  ela_pp->Print();
  TH1F * phi50 = new TH1F ("phi50","50 MeV Azimuthal angle",360,-180,180);
  TH1F * thetac50 = new TH1F ("thetac50","50 MeV CM Polar angle",180,0,180);
  ela_pp.Do(phi50,"p1 = [p,1]; p1->Boost([p + p]); _x = p1->Phi()* 180./TMath::Pi()");
  ela_pp.Do(thetac50,"p1 = [p,1]; p1->Boost([p + p]); _x = p1->Theta()* 180./TMath::Pi()");
  ela_pp.Loop(100000);
  
  test_elastic_phi->cd();
  phi50->SetMinimum(0);
  phi50->Draw("");
  test_elastic_theta->cd(1);
  thetac50->Draw("");

  PReaction ela_pp("0.030","p","p","p p");
  ela_pp->Print();
  TH1F * thetac30 = new TH1F ("thetac30","30 MeV CM Polar angle",180,0,180);
  ela_pp.Do(thetac30,"p1 = [p,1]; p1->Boost([p + p]); _x = p1->Theta()* 180./TMath::Pi()");
  ela_pp.Loop(100000);
  test_elastic_theta->cd(2);
  thetac30->Draw("");

  PReaction ela_pp("0.015","p","p","p p");
  ela_pp->Print();
  TH1F * thetac5 = new TH1F ("thetac15","15 MeV CM Polar angle",180,0,180);
  ela_pp.Do(thetac5,"p1 = [p,1]; p1->Boost([p + p]); _x = p1->Theta()* 180./TMath::Pi()");
  ela_pp.Loop(100000);
  test_elastic_theta->cd(3);
  thetac15->Draw("");

  PReaction ela_pp("0.005","p","p","p p");
  ela_pp->Print();
  TH1F * thetac5 = new TH1F ("thetac5","5 MeV CM Polar angle",180,0,180);
  ela_pp.Do(thetac5,"p1 = [p,1]; p1->Boost([p + p]); _x = p1->Theta()* 180./TMath::Pi()");
  ela_pp.Loop(100000);
  test_elastic_theta->cd(4);
  thetac5->Draw("");


}
