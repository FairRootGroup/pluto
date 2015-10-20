// Author: I. Froehlich
// Written: 1.02.2008
// Revised: 

#ifndef _PADAPTIVEMESHN_H_
#define _PADAPTIVEMESHN_H_

#define MAX_SUBMESH_N 20
#define MAX_DIMENSIONS 4

#include "TObject.h"
#include "PChannelModel.h"

class PAdaptiveMeshN : public TObject {

 private:

    Double_t y_max;        //upper value in each hyper-cube
    Double_t area_size;    //volume of the hyper-cube(es)
    Double_t x_max[MAX_DIMENSIONS], x_min[MAX_DIMENSIONS]; //width of the hypercubes
    Int_t sub_size[MAX_DIMENSIONS];       //sub-size in each dimension of the parent hypercube
    Int_t total_sub_size;  //total number of sub-meshs
    Int_t is_divided;
    Int_t layer;           //layer of the hypercube. The innermost has 0


    PAdaptiveMeshN *sub_tree;  //for each dimension pointer to the subtree
    Double_t *sub_area;        //for each dimension volume size

    Int_t num_dimensions;   //number of dimensions of the hyper-cube
    Int_t variable_dimensions;  //number of dimensions to be sampled
    PChannelModel * model;  //Model to be sampled
    Double_t array[MAX_DIMENSIONS]; //storage array for GetWeight()

    void ReCalc(void);      //re-calculate all volumes
    void FixArray(void);    //set the fixed array contributions
    

    UInt_t pattern;
    Int_t mcpoints;
    Double_t threshold_diff;
    Double_t threshold_abs;

    PAdaptiveMeshN * GetRandomBin(Double_t f_random);

 public:

    //constructor
    //pattern sets the variables to be sampled in the GetWeight method in PChannelModel
    PAdaptiveMeshN(UInt_t my_pattern, Int_t my_max_dimensions, 
		   PChannelModel * my_model, Double_t my_y_max);
    PAdaptiveMeshN();
    void SetDefaults(UInt_t my_pattern, Int_t my_max_dimensions, 
		     PChannelModel * my_model, Double_t my_y_max);
    void ReCalcYMax(void);
    ~PAdaptiveMeshN();

    void SetRange (Int_t dimension, Double_t my_x_min,  Double_t my_x_max) {
	x_min[dimension]=my_x_min;
	x_max[dimension]=my_x_max;
    }

    Bool_t GetRandom();

    Double_t GetArea() {return area_size;};
    
    Double_t GetArrayValue(Int_t dimension)  {
	if (dimension>num_dimensions) {
	    Warning("GetArrayValue","dimension>num_dimensions");
	    return 0;
	}
	return array[dimension];
    }

    Double_t GetXMax(Int_t dimension)  {
	// method should be fast. I do not make any test here!
	return x_max[dimension];
    };

    Double_t GetXMin(Int_t dimension) {
	// method should be fast. I do not make any test here!
	return x_min[dimension];
    };

    Double_t GetYMax() {
	return y_max;
    };


    void SetYMax(Double_t  d) {y_max=d;};

    void SetMCPoints(Int_t n) {mcpoints=n;};
    
    void SetThreshold(Double_t diff, Double_t abs) {
	threshold_diff=diff;
	threshold_abs=abs;
    };

    void Divide(Int_t num, Int_t my_layer);

    //void Draw(const Option_t *delme=NULL);

    void PrintMesh(void);

    ClassDef(PAdaptiveMeshN,0)  // Envelope with adaptable size of the bins (N dimensions)
};

#endif
