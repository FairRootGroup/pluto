void DrawAllFermiDistributions()
{
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetPalette(1);

   TStopwatch timer;
   timer.Start();
   double t0=-1;
   const char * percent="\%";
   int cpc1=1;
   int cpc2=1;
   int cpc3=1;
   int cpc4=1;
   int cpc5=1;
   int cpc6=1;

   //   const int Iterations = 100000000;
   const int Iterations = 1000000;
   const int Counter    = Iterations/100;

   double px1, py1, pz1;
   double px2, py2, pz2;
   double px3, py3, pz3;
   double px4, py4, pz4;
   double px5, py5, pz5;
   double px6, py6, pz6;

   TH1F * hD    = new TH1F("hD"   ,"Fermi Momentum Distribution inside D"   ,1000,0.,1.);
   TH1F * h3He  = new TH1F("h3He" ,"Fermi Momentum Distribution inside 3He" ,1000,0.,1.);
   TH1F * h4He  = new TH1F("h4He" ,"Fermi Momentum Distribution inside 4He" ,1000,0.,1.);
   TH1F * h7Li  = new TH1F("h7Li" ,"Fermi Momentum Distribution inside 7Li" ,1000,0.,1.);
   TH1F * h12C  = new TH1F("h12C" ,"Fermi Momentum Distribution inside 12C" ,1000,0.,1.);
   TH1F * h40Ca = new TH1F("h40Ca","Fermi Momentum Distribution inside 40Ca",1000,0.,1.);

   TCanvas *c1 = new TCanvas("c1","Fermi Distribution inside different Nuclei",200,10,600,400);

   TLegend *l1 = new TLegend(0.7,0.6,0.89,0.89);
   l1->SetFillColor(0);
   l1->AddEntry(hD,"D","l");
   l1->AddEntry(h3He,"3He","l");
   l1->AddEntry(h4He,"4He","l");
   l1->AddEntry(h7Li,"7Li","l");
   l1->AddEntry(h12C,"12C","l");
   l1->AddEntry(h40Ca,"40Ca","l");

   t0=timer.RealTime();
   
   makeDistributionManager()->Exec("nucleus_fermi:gamma");
   makeDistributionManager()->LinkDB();
   //------------------------------------------------------------------------------

   //    PFermiDistributions *fD    = new PFermiDistributions("d");
   
   //    for (int i=0; i<Iterations; i++)
   //    { 
   //       if (i%Counter==0)
   //       {
   //          printf(" %i%s of deuteron done in %f sec\n",cpc1-1,percent,timer.RealTime()-t0);
   //          ++cpc1;
   //          timer.Continue();
   //       }
   //       hD   ->Fill(fD   ->GetRandomFermiMomentum(px1,py1,pz1));
   //}
   //------------------------------------------------------------------------------

   PFermiMomentumGA* f3He = 
       ((PFermiMomentumGA*)makeDistributionManager()->GetDistribution("gp_in_3He"));   
   
   for (int i=0; i<Iterations; i++)
       { 
      if (i%Counter==0)
      {
         printf(" %i%s of helium-3 done in %f sec\n",cpc2-1,percent,timer.RealTime()-t0);
         ++cpc2;
         timer.Continue();
      }
      h3He ->Fill(f3He ->GetRandomFermiMomentum(px2,py2,pz2));
   }
   //------------------------------------------------------------------------------

   PFermiMomentumGA* f4He = 
       ((PFermiMomentumGA*)makeDistributionManager()->GetDistribution("gn_in_4He"));   
   for (int i=0; i<Iterations; i++)
   { 
      if (i%Counter==0)
      {
         printf(" %i%s of helium-4 done in %f sec\n",cpc3-1,percent,timer.RealTime()-t0);
         ++cpc3;
         timer.Continue();
      }
      h4He ->Fill(f4He ->GetRandomFermiMomentum(px3,py3,pz3));
   }
   //------------------------------------------------------------------------------
   
   PFermiMomentumGA* f7Li = 
       ((PFermiMomentumGA*)makeDistributionManager()->GetDistribution("gn_in_7Li"));   

   for (int i=0; i<Iterations; i++)
   { 
      if (i%Counter==0)
      {
         printf(" %i%s of lithium-7 done in %f sec\n",cpc4-1,percent,timer.RealTime()-t0);
         ++cpc4;
         timer.Continue();
      }
      h7Li ->Fill(f7Li ->GetRandomFermiMomentum(px4,py4,pz4));
   }
   //------------------------------------------------------------------------------
   
   PFermiMomentumGA* f12C = 
       ((PFermiMomentumGA*)makeDistributionManager()->GetDistribution("gn_in_12C"));   

   for (int i=0; i<Iterations; i++)
   { 
      if (i%Counter==0)
      {
         printf(" %i%s of carbon-12 done in %f sec\n",cpc5-1,percent,timer.RealTime()-t0);
         ++cpc5;
         timer.Continue();
      }
      h12C ->Fill(f12C->GetRandomFermiMomentum(px5,py5,pz5));
   }
   //------------------------------------------------------------------------------

   PFermiMomentumGA* f40Ca = 
       ((PFermiMomentumGA*)makeDistributionManager()->GetDistribution("gn_in_40Ca"));   

   for (int i=0; i<Iterations; i++)
   { 
      if (i%Counter==0)
      {
         printf(" %i%s of calcium-40 done in %f sec\n",cpc6-1,percent,timer.RealTime()-t0);
         ++cpc6;
         timer.Continue();
      }
      h40Ca->Fill(f40Ca->GetRandomFermiMomentum(px6,py6,pz6));
   }
   //------------------------------------------------------------------------------

   //   hD->Scale(1./hD->GetMaximum());
   h3He->Scale(1./h3He->GetMaximum());
   h4He->Scale(1./h4He->GetMaximum());
   h7Li->Scale(1./h7Li->GetMaximum());
   h12C->Scale(1./h12C->GetMaximum());
   h40Ca->Scale(1./h40Ca->GetMaximum());


   c1->cd(1);
   hD   ->SetLineColor(1);
   hD   ->SetLineWidth(2);
   h3He ->SetLineColor(2);
   h3He ->SetLineWidth(2);
   h4He ->SetLineColor(3);
   h4He ->SetLineWidth(2);
   h7Li ->SetLineColor(4);
   h7Li ->SetLineWidth(2);
   h12C ->SetLineColor(5);
   h12C ->SetLineWidth(2);
   h40Ca->SetLineColor(6);
   h40Ca->SetLineWidth(2);


   hD->SetStats(kFALSE);
   hD->SetLabelSize(.04,"X");
   hD->SetLabelSize(.04,"Y");
   hD->SetTitleSize(.04,"X");
   hD->SetTitleSize(.04,"Y");
   hD->SetLabelFont(41,"X");
   hD->SetLabelFont(41,"Y");
   hD->SetNdivisions(505,"X");
   hD->SetNdivisions(505,"Y");
   hD->GetXaxis()->CenterTitle(kTRUE);
   hD->SetXTitle("#font[41]{ Fermi Momentum p_{f} [GeV]}");
   hD->GetYaxis()->CenterTitle(kTRUE);
   hD->SetYTitle("#font[41]{Counts [a.u.]}");
   hD->SetTitleOffset(1.3,"X");
   hD->SetTitleOffset(1.,"Y");
   hD->Draw();
   h3He ->Draw("same");
   h4He ->Draw("same");
   h7Li ->Draw("same");
   h12C ->Draw("same");
   h40Ca->Draw("same");

   l1->Draw();
   c1->Update();

}
