// Author: Romain Holzmann
// Written: 8/11/2002
// Modified: 8/11/2002 R.H.
// PFileInput Class Header

#ifndef _PFileInput_H_
#define _PFileInput_H_

#include "PParticle.h"
#include "PChannel.h"

class PFileInput: public PParticle {

private:

 protected:
  Int_t imode;      //  mode (HGeant=1, QMD=2, BUU=3)
  Char_t* name;     //! input generator name
  Char_t* file;     //! input file name
  PParticle** part; //! array of particle pointers
  FILE* fp;         //! input file pointer
  Int_t flag;       //! event header flag

 public:
  PFileInput(Char_t* mode, Char_t* filename);

  void setToMidrapidity(float agev);
  Int_t isFileInput() { return 1; } 
  Int_t readEventHeader(Float_t& b);
  Int_t readParticle(Double_t& px, Double_t& py, Double_t& pz, Double_t& E,
                     Double_t& vx, Double_t& vy, Double_t& vz, Double_t& vt,
                     Int_t& Id, Int_t& srcId, Int_t& parId, Int_t& parInd, Double_t& weight);
  PChannel* makeChannel(Int_t nMax, Float_t Ebeam=0.);
  virtual void Print(const Option_t* delme=NULL) const;
  virtual ~PFileInput() {
    if (fp) fclose(fp);
  }

ClassDef(PFileInput,1) // Pluto file input
};
#endif // _PFileInput_H_

















