///////////////////////////////////////////////////////////////////////////////
//  Pluto Utilities Class
//
//  This is a static class containing utility functions.
//
//                             Author:  M.A. Kagarlis
//                             Written: 15.12.98
//                             Revised: 15.06.00
//
// Ref 1: CERNLIB function FLPSOR, R.S. Scowen, Algorithm 271 QUICKERSORT,
//        Collected Algorithms from CACM (1965) (sorting function dsort)
// Ref 2: D.M. Brink and G.R. Satchler, Angular Momentum, Oxford Library (1961)
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TMath.h"


const Double_t fact[34]={1.,1.,2.,6.,2.4e1,1.2e2,7.2e2,5.04e3,4.032e4,
		       3.6288e5,3.6288e6,3.99168e7,4.790016e8,6.2270208e9,
		       8.71782912e10,1.307674368e12,2.0922789888e13,
		       3.556874281e14,6.4023737057e15,1.2164510041e17,
		       2.4329020082e18,5.1090942172e19,1.1240007278e21,
		       2.5852016739e22,6.2044840173e23,1.5511210043e25,
		       4.0329146113e26,1.088886945e28,3.0488834461e29,
		       8.8417619937e30,2.6525285981e32,8.2228386542e33,
		       2.6313083693e35,8.6833176188e36}; // factorials to 33!

#include "PUtils.h"


void PUtils::correct_histo(TH1 * histo) {
    Int_t nbins=0;
    for(int i=0;i<histo->GetNbinsX() ; i++)
	if (histo->GetBinContent(i)) nbins++;
    histo->Scale(nbins);
}

void PUtils::correct(TH1 * histo) {
    histo->Scale(1./histo->GetBinWidth(1));
}

void PUtils:: dsort(Double_t * a, int n) {
// Sort in ascending order the first n entries of the array *a.
  int i,j,l,m,r,rt[20],lt[20],level=1;
  Double_t x,w;
  lt[0]=1;
  rt[0]=n;
 a: 
  l=lt[level-1];
  r=rt[level-1];
  --level;
 b:
  if (r<=l) {
    if (level>0) goto a;
    else return;
  }
  i=l;
  j=r;
  m=(l+r)/2;
  x=*(a+m-1);
//  cout << "X: " << x << endl;
 c:
  while (*(a+i-1)<x) ++i;
  while (*(a+j-1)>x) --j;
  while (i<=j) {
    w=*(a+i-1);
    *(a+i-1)=*(a+j-1);
    *(a+j-1)=w;
    ++i;
    --j;
    if (i<j) goto c;
  }
  ++level;
  if (r-i>=j-l) {
    lt[level-1]=i;
    rt[level-1]=r;
    r=j;
  } else {
    lt[level-1]=l;
    rt[level-1]=j;
  }
  goto b;
}

Double_t PUtils:: cgc(const int & jx1, const int & jx2, const int & jx3,
		    const int & mx1, const int & mx2) {
  // Clebsch-Gordan coefficients (arguments are 2x the corresponding j and m)
  int j1=jx1, j2=jx2, j3=jx3, m1=mx1, m2=mx2, icntr, it, jz1, jz2, jz3,
    jt1, jt2, jt3, jt4, jt5, j4, j5, numin, numax;
  Double_t vcc=0., fctor, phas;
  if ((j1+j2+j3)%2||(j1+m1)%2||(j2+m2)%2) return vcc;
  j1=jx1;
  j2=jx2;
  j3=jx3;
  m1=mx1;
  m2=mx2;
  if (j1<j2) goto a;
  if (j3<j2) goto b;
  icntr=0;
  goto c;
 a:
  if (j3<j1) goto b;
  icntr=-1;
  it=j1;
  j1=j2;
  j2=it;
  it=m1;
  m1=m2;
  m2=it;
  goto c;
 b:
  icntr=1;
  it=j2;
  j2=j3;
  j3=it;
  m2=-m1-m2;
 c:
  jz1=(j1+j2-j3)/2;
  if (jz1<0) return vcc;
  jz2=(j1+j3-j2)/2;
  if (jz2<0) return vcc;
  jz3=(j2+j3-j1)/2;
  if (jz3<0||j1-abs(m1)<0||j2-abs(m2)<0||j3-abs(m1+m2)<0) return vcc;
  jt1=(j1-j3+m2)/2;
  jt2=(j2-j3-m1)/2;
  numin=TMath::Max(TMath::Max(jt1,jt2),0);
  jt3=(j1-m1)/2;
  jt4=(j2+m2)/2;
  numax=TMath::Min(TMath::Min(jt3,jt4),jz1);
  jt5=(j2-m2)/2;
  if (numax<numin) return vcc;
  j4=j1/2;
  j5=j3/2;
  phas=phasef(numin);
  for (int nu=numin;nu<=numax;++nu) {
    vcc+=phas*fact[j4]*fact[j5]/
      (fact[jt3-nu]*fact[nu-jt2]*fact[jt4-nu]*
       fact[nu-jt1]*fact[jz1-nu]*fact[nu]);
    phas=-phas;
  }
  fctor=fact[(j1+m1)/2]*fact[jt3]*fact[jz2]*fact[(j3+m1+m2)/2]*
    fact[(j3-m1-m2)/2]*fact[jz1]*fact[jz3]*fact[jt4]*fact[jt5]*(j3+1)/
    (fact[j4]*fact[j4]*fact[(j1+j2+j3)/2+1]*fact[j5]*fact[j5]);
  vcc*=sqrt(fctor);
  if (icntr<0) vcc*=phasef(jz1);
  else if (icntr>0) vcc*=sqrt(Double_t(j2+1)/Double_t(j3+1))*phasef(jt3);
  return vcc;
}

Double_t PUtils:: s3j(const Double_t & jx1, const Double_t & jx2, const Double_t & jx3,
		    const Double_t & jm1, const Double_t & jm2, Double_t jm3) {
  // 3j-symbol, related to Clebsch-Gordan coefficient
  if (jm3!=0&&jm1+jm2+jm3!=0.) return 0.;
  int j1=2*int(jx1), j2=2*int(jx2), j3=2*int(jx3),
    m1=2*int(jm1), m2=2*int(jm2);
  return cgc(j1,j2,j3,m1,m2)*phasef((j1-j2+m1+m2)/2)/sqrt(1.+j3);
}

Double_t PUtils:: racah(const int & j1, const int & j2, const int & j3,
		      const int & j4, const int & j5, const int & j6) {
  // Racah coefficients (arguments are 2 x j or m)
  Double_t z1=j123(j1,j2,j5), z2=j123(j1,j3,j6);
  if (z1==0.||z2==0.) return 0.;
  z1*=j123(j3,j4,j5);
  if (z1==0.) return 0.;
  z2*=j123(j2,j4,j6);
  if (z2==0.) return 0.;
  z1=z2*sqrt(z1/z2);
  int jt1=(j1+j2+j5)/2, jt2=(j3+j4+j5)/2, jt3=(j1+j3+j6)/2,
    jt4=(j2+j4+j6)/2, jz1=(j1+j2+j3+j4)/2, jz2=(j1+j4+j5+j6)/2,
    jz3=(j2+j3+j5+j6)/2, numin=TMath::Max(TMath::Max(jt1,jt2),TMath::Max(jt3,jt4)),
    numax=TMath::Min(TMath::Min(jz1,jz2),jz3);
  if (numax<numin) return 0.;
  Double_t phase=phasef(numin+jz1)*z1, r=0.;
  for (int nu=numin;nu<=numax;++nu) {
    int jy1=nu-jt1, jy2=nu-jt2, jy3=nu-jt3,
      jy4=jz1-nu, jy5=jz2-nu, jy6=jz3-nu;
    Double_t fctor=fact[jy1]*fact[jy2]*fact[jy3]*fact[nu-jt4]*
      fact[jy4]*fact[jy5]*fact[jy6]/fact[nu+1];
    if (fctor>0.) r+=phase/fctor;
  }
  return r;
}

Double_t PUtils:: s6j(const Double_t & x1, const Double_t & x2, const Double_t & x3,
		    const Double_t & x4, const Double_t & x5, const Double_t & x6) {
  // 6j-symbol, related to Racah coefficient
  int j1=2*int(x1), j2=2*int(x2), j3=2*int(x3),
    j4=2*int(x4), j5=2*int(x5), j6=2*int(x6);
  return racah(j1,j2,j5,j4,j3,j6)*phasef((j1+j2+j4+j5)/2);
}

Double_t PUtils:: j123(const int & j1, const int & j2, const int & j3) {
  // used by Racah
  int jz1=(j1+j2-j3)/2, jz2=(j1-j2+j3)/2, jz3=(j2+j3-j1)/2, jz4=(j1+j2+j3)/2+1;
  return (jz1<0||jz2<0||jz3<0) ? 0. : fact[jz1]*fact[jz2]*fact[jz3]/fact[jz4];
}
  
Int_t PUtils::FindIndex(Int_t n, Double_t *a, Double_t r) {
  // return index i of first element a[i] which is larger than r.
  // a[] is assumed to be sorted in ascending order
  Int_t imin=0;
  Int_t imax = n-1;
  Int_t ihalf = (imin+imax)/2;
  if (a[0]>r) return 0;
  while (ihalf>imin) {  // algorithm is by binary search
    if (a[ihalf]<r) imin=ihalf;
    else imax=ihalf;
    ihalf=(imin+imax)/2;
  }
  return imax;

}

//TRandom3 PUtils::REngine;

Bool_t PUtils::Tokenize(const char * options, const char *delimiter, char ** array, int *size) {
    //Tokenize the string "options" using
    //the "delimiter" 
    //(like "+" or ",")
    //The char array is constructed but never deleted
    //in addition, any spaces before and after the delimiter are cleaned
    //-> for setup phase, but not for event loops
    //The user has to provide the pointer to the char array pointers
    //with "size", these pointers are updated
    //In addition, the size is also updated
    //
    //The ROOT version is buggy:
    //TString spattern("a - b -> c");
    //TObjArray *carray = spattern.Tokenize(TString("->"));     
    //carray->GetEntriesFast()
    //(const Int_t)3
    

    //cout << options << " -> " << delimiter << endl;
    for (int i=0;i<*size;i++) array[i]=NULL;
    
    char * mystack=new char[strlen(options)+1];
    strcpy(mystack,options);

    Int_t pat=1;
    array[0]=mystack;

    while (strstr(mystack,delimiter)) {
	if (pat==*size) {
	    cout << "PUtils::Tokenize: size " << *size<< " is too small for " << options << " ,delim:" << delimiter << endl;
	}
	char * pos = strstr(mystack,delimiter);
	*pos='\0';
	mystack = pos + strlen(delimiter);
	array[pat]=mystack;
	//cout << "+-: "<<array[pat] << endl;
	pat++;
    }

    *size=pat;

    for (int i=0;i<*size;i++) {
	//remove leading/trailing spaces

	remove_spaces(&array[i]);

    }

    return kTRUE;

}


void PUtils::remove_spaces(char ** partc) {
    
    while (**partc==' ') (*partc)++;
    if (strlen(*partc)) {
	int partend=strlen(*partc)-1;
//	    cout << partc << ":" << partc[partend] << endl;
	while ((*partc)[partend]==' ' && partend>=0) {
	    (*partc)[partend]='\0';
	    partend--;
	}
    }
}


Int_t PUtils::remove_brackets(char ** partc, char a, char b) {
    //remove matching brackets like a=(, b=)
    //cout << "in: " << *partc << endl;
    remove_spaces(partc);

    if ((*partc)[0] != a) return 0;

    Int_t loop=1;
    Int_t found=0;
    Int_t max_br=0;

    while (loop) {
	loop=0;

	if (((*partc)[0] == a) && ((*partc)[(strlen(*partc)-1)] == b)) {
	    //brackets at begin and end

	    //do they match?
	    Int_t match=1;
	    Int_t brackets=0;
	    for (UInt_t i=0;i<strlen(*partc);i++) {

		if ((*partc)[i] == a) { 
		    max_br++;
		    brackets++;
		}
		if (((*partc)[i] == b) && max_br ) {
		    brackets--;
		    if (!brackets && i != (strlen(*partc)-1)) {
			//matching bracket is not last character
			match=0;
		    }
		}
	    }
	    if (match) {
		found++;
		(*partc)[(strlen(*partc)-1)]='\0';
		(*partc)++;
		
		remove_spaces(partc);
		loop=1;
	    }
	    
	}

    }




    //cout << "out: " << *partc << endl;

    return found;

}

Bool_t PUtils::ValidVariableName(const char * name) {
    //We check if "name" is a variable name
    //or a possible command
    //Allowed names for variables is alphanumeric + "_"
    //First char must be alpha or "#"

    Bool_t isvar = kTRUE;
    if (!(name[0]=='#') and !(name[0]=='_') and !isalpha(name[0])) isvar = kFALSE;

    for (unsigned int i=1; i<strlen(name); i++)
	if (!(name[i]=='_') and !isalnum(name[i])) isvar = kFALSE;

    return isvar;
}


Bool_t PUtils::IsInt(const char * name) {
    //We check if "name" is an integer variable name

    Bool_t isvar = kTRUE;

    for (unsigned int i=1; i<strlen(name); i++)
	if (!isdigit(name[i])) isvar = kFALSE;
    
    return isvar;
}

PUtilsREngine& fPUtilsREngine() {
    static PUtilsREngine* ans = new PUtilsREngine();
    return *ans;
}

PUtilsREngine* makePUtilsREngine() {
    return &fPUtilsREngine();
}


PUtilsREngine::PUtilsREngine () {  
    rnd = new TRandom3();     
    SetSeed(SEED);  
} 

void PUtilsREngine::SetSeed(UInt_t s) { 
    rnd->SetSeed(s); 

    if (s)
	Warning("PUtilsREngine","Seed set FIXED to %i",s);
    else
	Info("PUtilsREngine","Random seed set to %i",rnd->GetSeed());
}

ClassImp(PUtils)
ClassImp(PUtilsREngine)
