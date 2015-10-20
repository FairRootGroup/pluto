////////////////////////////////////////////////////////
//  An adaptive mesh for enveloping the TF1 function
//  Used for fast random sapmling when the ROOT 
//  GetRandom function cannot be used
//  The build-in ROOT random sampling has really big
//  disadvantages when parameters (e.g. kinetic energy of the 
//  parent) are not stable and the shape of the function is changes
//  In this case ROOT re-creates the polymonials and this 
//  takes AGES!
//
//  The adaptive mesh method here is based on the rejection
//  meshod, the test function is a step function which is 
//  adaptable: If the function is steep, more bins are created
//  
//
//                    Author: I. Froehlich
//                    Written: 1.02.2003
//                    Revised: 
//
////////////////////////////////////////////////////////


#include "PAdaptiveMesh.h"
#include "PUtils.h"
#include <iostream>
#include "TH1.h"
#include <cmath>


using namespace std;

PAdaptiveMesh::PAdaptiveMesh(Double_t my_x_min, Double_t my_x_max, Double_t my_y_max) {
    x_max=my_x_max;
    x_min=my_x_min;
    y_max=my_y_max;
    sub_size=layer=0;
    tf1=NULL;
    line=NULL;

    area_size=(x_max-x_min)*y_max;


}

PAdaptiveMesh::~PAdaptiveMesh() {

};

void PAdaptiveMesh::ReCalc() {
    
    if (sub_size) {
	area_size=0;
	for (int i=0;i<sub_size;i++) {
	    sub_tree[i]->ReCalc();
	    area_size+=sub_tree[i]->GetArea();
	    sub_area[i]=area_size;

	    //   cout << layer << ":"<<sub_tree[i]->GetArea() << ":" << sub_area[i] << endl;
	}
    } else {
	area_size=(x_max-x_min)*y_max;
	//cout << x_min << "," << x_max<< endl;
    }

    

};


PAdaptiveMesh * PAdaptiveMesh::GetRandomBin(Double_t f_random) {
    if (!sub_size) return this;
    
    for (int i=0;i<sub_size;i++) {
//	    cout << total_sub_size << endl;
	if (f_random < sub_area[i]) {
//		cout << "b" << endl;
	    if (i) return sub_tree[i]->GetRandomBin(f_random - sub_area[i-1]);
	    return sub_tree[i]->GetRandomBin(f_random);

	}
    }
    Fatal("GetRandomBin","Reached and of function");

    return NULL;

}

Double_t PAdaptiveMesh::GetRandom() {
    //First check which bin we have to access


 repeat2:

    PAdaptiveMesh * bin=GetRandomBin(PUtils::sampleFlat() * area_size);

    Double_t x_random=PUtils::sampleFlat()*(
	bin->GetXMax()-bin->GetXMin()
	) + bin->GetXMin();
    
    
    Double_t y_random=tf1->Eval(x_random);

    // cout << x_random << ":" << y_random << endl;

    if (y_random>bin->GetYMax()) {
	//bin has a maximum larger then expected
	//re-scale everything!
	bin->SetYMax(y_random);

	cout << "recalc" << endl;

	ReCalc();
	goto repeat2;

    }

    if (y_random > (bin->GetYMax()*PUtils::sampleFlat())) {
// cout << x_random <<endl;
	return x_random;
    }
    
//    cout << x_random << ":" << bin->GetYMax() << ":" << bin->GetArea()  << endl;

//    exit(1);

    goto repeat2;

}


void PAdaptiveMesh::Divide(Int_t num, Int_t my_layer) {
    if (sub_size) return; //already divided

    if (num > MAX_SUBMESH) {
	Warning("Divide","num > MAX_SUBMESH");
	return;
    }

    if (!tf1) {
	Warning("Divide","No TF1 function set");
	return;
    }

    sub_size=num;
    for (int i=0;i<sub_size;i++) {
	Double_t min=x_min+((x_max- x_min)/(Double_t)num)*(Double_t)i;
	Double_t max=min+((x_max- x_min)/(Double_t)num);
	Double_t ymax=tf1->Eval(min);
	if (tf1->Eval(max)> ymax) ymax=tf1->Eval(max);

	if (ymax<0) ymax=0;

//	cout << "-->" << min << ":" << max << endl;

	sub_tree[i]=new PAdaptiveMesh(min,max,ymax);
	sub_tree[i]->SetTF1(tf1);
	
	layer=my_layer;

	if (layer>0) {

	    Double_t diff = tf1->Eval(min)/tf1->Eval(max);
	    Double_t abs_diff = fabs(tf1->Eval(min)-tf1->Eval(max));

	    if (((diff > 1.1) || ((1./diff)>1.1)) &&
		(abs_diff>0.5)) {
		sub_tree[i]->Divide(num,layer-1);
	    }

	}
	


    }

    ReCalc();

}

void PAdaptiveMesh::Draw(const Option_t*) {
    if (sub_size) {
	for (int i=0;i<sub_size;i++) {
	    sub_tree[i]->Draw();
	}
	
    } else {
	if (!line)
	    line = new TLine(x_min, y_max, x_max, y_max);
	else {
	    line->SetX1(x_min);
	    line->SetX2(x_max);
	    line->SetY1(x_max);
	    line->SetY2(x_max);
	}
	line->Draw();
    }
    tf1->Draw();
}
