// plot pt vs. y lines on top of a graph
//
void plot(Float_t m, Float_t y0 = 0) {
  for (Float_t p=0.1; p<2.0; p+=0.1) {
    Float_t ymax = TMath::ASinH(p/m);  // compute extrema
    Float_t ymin = -ymax;
    ymin += y0;                        // shift from cm to lab frame
    ymax += y0;
    TF1 *f = new TF1("pty","1000*sqrt(([0]**2+[1]**2)/cosh(x-[2])**2-[0]**2)",ymin,ymax); 
    f->SetParameters(m,p,y0);
    f->SetLineWidth(1);
    f->SetLineColor(2);
    f->SetLineStyle(2);
    f->DrawCopy("same");
    delete f;
  }
}
    f->SetLineColor(2);
