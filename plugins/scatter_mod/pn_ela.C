//TITLE Low-energy pp scattering

{
  gROOT->Reset();

  makeDistributionManager()->Enable("pn_elastic");
  
  TCanvas *test_elastic_phi = new
    TCanvas ("test_elastic_phi","pn elastic",100,10,600,400);
  
  TCanvas *test_elastic_theta = new
    TCanvas ("test_elastic_theta","pn elastic",200,5,700,800);
  test_elastic_theta->Divide(2,2);
  
  PReaction ela_pn("0.0011","p","n","p n");
  ela_pn->Print();
  TH1F * phi1 = new TH1F ("phi10","1 MeV Azimuthal angle",360,-180,180);
  TH1F * thetac1 = new TH1F ("thetac1","1 MeV CM Polar angle",180,0,180);
  ela_pn.Do(phi1,"p1 = [p,1]; p1->Boost([p + n]); _x = p1->Phi()* 180./TMath::Pi()");
  ela_pn.Do(thetac1,"p1 = [p,1]; p1->Boost([p + n]); _x = p1->Theta()* 180./TMath::Pi()");
  ela_pn.Loop(100000);
  
  test_elastic_phi->cd();
  phi1->SetMinimum(0);
  phi1->Draw("");
  test_elastic_theta->cd(1);
  thetac1->Draw("");

  return;

  PReaction ela_pn("0.159","p","n","p n");
  ela_pn->Print();
  TH1F * thetac159 = new TH1F ("thetac159","159 MeV CM Polar angle",180,0,180);
  ela_pn.Do(thetac159,"p1 = [p,1]; p1->Boost([p + n]); _x = p1->Theta()* 180./TMath::Pi()");
  ela_pn.Loop(100000);
  test_elastic_theta->cd(2);
  thetac159->Draw("");

  PReaction ela_pn("0.370","p","n","p n");
  ela_pn->Print();
  TH1F * thetac370 = new TH1F ("thetac370","370 MeV CM Polar angle",180,0,180);
  ela_pn.Do(thetac370,"p1 = [p,1]; p1->Boost([p + n]); _x = p1->Theta()* 180./TMath::Pi()");
  ela_pn.Loop(100000);
  test_elastic_theta->cd(3);
  thetac370->Draw("");

  PReaction ela_pn("2.0","p","n","p n");
  ela_pn->Print();
  TH1F * thetac2000 = new TH1F ("thetac2000","2000 MeV CM Polar angle",180,0,180);
  ela_pn.Do(thetac2000,"p1 = [p,1]; p1->Boost([p + n]); _x = p1->Theta()* 180./TMath::Pi()");
  ela_pn.Loop(100000);
  test_elastic_theta->cd(4);
  thetac2000->Draw("");


}
