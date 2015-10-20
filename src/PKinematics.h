// Author: M.A. Kagarlis
// Written: 6.1.00
// PKinematics Class Header

#ifndef _PKINEMATICS_H_
#define _PKINEMATICS_H_

#include "PStaticData.h"

class PKinematics {
  // useful static kinematics functions

 public:
  static double lambda(double M, double m1, double m2) {
    // Kaellen function
    double y=(m1+m2);
    if (M<y) return 0.;
    y*=y;
    double x=M*M, z=(m1-m2)*(m1-m2), w=(x-y)*(x-z);
    return (w>0.) ? w : 0.;
  }

  static double pcms2(double M, double m1, double m2) {
    // cm momentum^2 for the decay of M to m1 and m2
    return (M>0.) ? 0.25*lambda(M,m1,m2)/(M*M) : 0.;
  }
  
  static double pcms(double M, double m1, double m2) {
    // cm momentum for the decay of M to m1 and m2
      return sqrt(pcms2(M,m1,m2));
  }
  
  static double pcmt(double m, double t) {
    const double mp=makeStaticData()->GetParticleMass("p"), 
	mp2=mp*mp;
    double m2=m*m;
    double arg=m2+mp2-t;
    arg=0.25*arg*arg/m2 - mp2;
    return (arg>0.) ? sqrt(arg) : 0. ;
  }

  static void rotes(double c, double s, double c2, double s2, double *pr) {
    // rotation used by GENBOD
    double sa=*pr, sb=*(pr+1), a=sa*c-sb*s;
    *(pr+1)=sa*s+sb*c;
    double b=*(pr+2);
    *pr=a*c2-b*s2;
    *(pr+2)=a*s2+b*c2;
  }
};
#endif // _KINEMATICS_H_
