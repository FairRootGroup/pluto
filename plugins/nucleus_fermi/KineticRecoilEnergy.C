
TLorentzVector doCalEnergy(double BeamEnergy,
                           TLorentzVector Particle1,
                           TLorentzVector Particle2,
                           double nucleusMass,
                           double Particle2Mass,
                           double Particle3Mass)
{
    double E_Particle1   = Particle1.E();
    double p_Particle1_x = Particle1.Px();
    double p_Particle1_y = Particle1.Py();
    double p_Particle1_z = Particle1.Pz();
    double p_Particle1   = sqrt(TMath::Power(p_Particle1_x,2.0) +
                                    TMath::Power(p_Particle1_y,2.0) +
                                    TMath::Power(p_Particle1_z,2.0));
    double phi   = Particle2.Phi();
    double theta = Particle2.Theta();
    double b     = 2.0 * ( p_Particle1_x * cos(phi) * sin(theta) +
                           p_Particle1_y * sin(phi) * sin(theta) +
                           p_Particle1_z * cos(theta) -
                           BeamEnergy * cos(theta)
                         );
    double c     = p_Particle1 * p_Particle1 + BeamEnergy * BeamEnergy - 2.0 * BeamEnergy * p_Particle1_z;
    double d     = BeamEnergy + nucleusMass - E_Particle1;
    double e     = TMath::Power(Particle3Mass,2.0) - TMath::Power(Particle2Mass,2.0) - d * d + c;
    double Delta = 16.0 * TMath::Power(d,2.0) * (TMath::Power(e,2.0) +
                                                 TMath::Power(b * Particle2Mass,2.0) -
                                                 TMath::Power(d * Particle2Mass * 2.0,2.0));
    
    TLorentzVector NewParticle(0.0,0.0,0.0,0.0);
    if(Delta>0.)
    {
       double sol2     = (2.0 * e * b + sqrt(Delta)) / (2.0 * (4.0 * TMath::Power(d,2.0) - TMath::Power(b,2.0)));
       double newpxcal = sol2 * cos(phi) * sin(theta);
       double newpycal = sol2 * sin(phi) * sin(theta);
       double newpzcal = sol2 * cos(theta);
       double energy   = sqrt(TMath::Power(sol2,2.0) + TMath::Power(Particle2Mass,2.0));
      
       TLorentzVector NewParticle2(newpxcal,newpycal,newpzcal,energy);
       NewParticle = NewParticle2;
    }
  
    return NewParticle;
}
