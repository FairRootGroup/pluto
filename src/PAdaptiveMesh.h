// Author: I. Froehlich
// Written: 1.02.2003
// Revised: 

#ifndef _PADAPTIVEMESH_H_
#define _PADAPTIVEMESH_H_

#define MAX_SUBMESH 20

#include "TObject.h"
#include "TF1.h"
#include "TLine.h"

class PAdaptiveMesh : public TObject {

 private:

    Double_t y_max;        //upper boundary
    Double_t area_size;    //area of the box(es)
    Double_t x_max, x_min; //bin width
    Int_t sub_size, layer;
    PAdaptiveMesh * sub_tree [MAX_SUBMESH];
    Double_t  sub_area [MAX_SUBMESH];
    TF1 * tf1;
    TLine *line;  //upper line for the Draw() function 
    void ReCalc(void);
    PAdaptiveMesh * GetRandomBin(Double_t f_random);
 public:

    //constructor
    PAdaptiveMesh( Double_t my_x_max, Double_t my_x_min, Double_t my_y_max);
    ~PAdaptiveMesh();
 
    void SetTF1(TF1 * t){tf1=t;};
    Double_t GetRandom();
    Double_t GetArea() {return area_size;};
    

    Double_t GetXMax() {return x_max;};
    Double_t GetXMin() {return x_min;};
    Double_t GetYMax() {return y_max;};
    void SetYMax(Double_t  d) {y_max=d;};

    void Divide(Int_t num, Int_t my_layer);

    void Draw(const Option_t *delme=NULL);

    ClassDef(PAdaptiveMesh,0)  //TF1 envelope with adaptable size of the bins
};

#endif
