////////////////////////////////////////////////////////
//  An adaptive mesh for enveloping the GetWeight() function
//  of a PChannelModel
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
//                    Written: 1.02.2008
//                    Revised: 
//
////////////////////////////////////////////////////////


#include "PAdaptiveMeshN.h"

#include "PUtils.h"
#include <iostream>
#include <cmath>

#include "TStopwatch.h"

using namespace std;


PAdaptiveMeshN::PAdaptiveMeshN(UInt_t my_pattern, Int_t my_max_dimensions, 
			       PChannelModel * my_model, Double_t my_y_max) {
    // Constructor for the adaptive mesh
    // The "pattern" defines which dimensions of the GetWeight method
    // of the PChannelModel "my_model" will be sampled. The LSB
    // coressponds to the lowest value in the array "mass"
    //
    // Fixed dimension should be set with SetXMax to the nomal value
    // The SetYMax value will eb ignored in this case
    // "my_y_max" is the starting value for the rejection limit.
    // It will be dynamically updated - but it should not be too far from the
    // max value of the GetWeight
    
    SetDefaults(my_pattern, my_max_dimensions, 
		my_model, my_y_max);
    
}

PAdaptiveMeshN::PAdaptiveMeshN() {
    
}

void PAdaptiveMeshN::SetDefaults(UInt_t my_pattern, Int_t my_max_dimensions, 
				 PChannelModel * my_model, Double_t my_y_max) {

    if (my_max_dimensions>MAX_DIMENSIONS) {
	Fatal("SetDefaults","my_max_dimensions>MAX_DIMENSIONS");
    }
    
    model = my_model;
    pattern = my_pattern;
    y_max=my_y_max;


    num_dimensions = my_max_dimensions;
    variable_dimensions = 0;

    for (int i=0;i<num_dimensions;i++) {
	if (((pattern >> i) & 0x1) == 0x1) { 
	    variable_dimensions++;
	    x_max[i]=0.;
	    x_min[i]=0.;
	} else {
	    x_min[i]=x_max[i]=0.;
	}
	sub_size[i]=0;
    }
    is_divided=0; 
    total_sub_size=0;
    threshold_abs=0.1;
    threshold_diff=1.1;
    mcpoints=0;

    sub_tree = NULL;
    sub_area = NULL;
    
    ReCalc();

};

PAdaptiveMeshN::~PAdaptiveMeshN() {
    
    // cout << "delete mesh" << endl;

};

void PAdaptiveMeshN::ReCalc() {
    
    area_size=0.;
    for (int i=0;i<total_sub_size;i++) {    

	sub_tree[i].ReCalc();
	area_size+=sub_tree[i].GetArea();
	sub_area[i]=area_size;	
	//cout << layer << ":"<<sub_tree[i].GetArea() << ":" << sub_area[i] << endl;
    }

    if (!area_size) {
	area_size=y_max;
	for (int i=0;i<num_dimensions;i++) {    
	    if (((pattern >> i) & 0x1) == 0x1) { 
		area_size*=(x_max[i]-x_min[i]);
	    }
	    //cout << x_min << "," << x_max<< endl;
	}
    }

};


PAdaptiveMeshN * PAdaptiveMeshN::GetRandomBin(Double_t f_random) {
    if (!total_sub_size) return this;
    
    for (int i=0;i<total_sub_size;i++) {
	//	    cout << total_sub_size << endl;
	if (f_random < sub_area[i]) {
//		cout << "b" << endl;
	    if (i) return sub_tree[i].GetRandomBin(f_random - sub_area[i-1]);
	    return sub_tree[i].GetRandomBin(f_random);

	}
    }
    Error("GetRandomBin","Reached end of function");
    PrintMesh();
    
    return NULL;

}

Bool_t PAdaptiveMeshN::GetRandom() {
    //First check which bin we have to access
    Int_t num=0;

    TStopwatch timer;                        // time loop
    timer.Start();
 
 repeat2:
    //num=0;
    //cout << "a" << endl;
    PAdaptiveMeshN * bin=GetRandomBin(PUtils::sampleFlat()*area_size);
    //PAdaptiveMeshN * bin=GetRandomBin(0.);
    //cout << "b" << endl;
 repeat3:
    num++;

    //Fill a random number somewhere in the sub-bin
    FixArray();
    for (int i=0;i<num_dimensions;i++) {
	if (((pattern >> i) & 0x1) == 0x1) { 
	    array[i]=
		PUtils::sampleFlat()*(bin->GetXMax(i)-bin->GetXMin(i)) 
		+ bin->GetXMin(i);
	}
    }

    return kTRUE;

    // cout << array[0] << ":" << array[1] << endl;
    //Double_t y_random = model->GetWeight(array);
    Double_t y_random = model->GetWeight(array);
    //Double_t y_random = 0;
    //y_random = 0;

    //cout << y_random << "->" << bin->GetYMax() << endl;

    if (y_random>bin->GetYMax()) {
 	//bin has a maximum larger then expected
 	//re-scale everything!
 	bin->SetYMax(y_random);

 	cout << "recalc" << endl;

 	ReCalc();
 	goto repeat2;

    }

     
    //if (num>10) return kTRUE;

    if (y_random > (bin->GetYMax()*PUtils::sampleFlat())) {
// // cout << x_random <<endl;
//	cout << "num: " << num << endl;
 	return kTRUE;
    }
    
//    cout << x_random << ":" << bin->GetYMax() << ":" << bin->GetArea()  << endl;

//    exit(1);

   
    //printf("%f sec\n",timer.RealTime());timer.Continue();

    if (num<100) goto repeat3;
    goto repeat2;

}


void PAdaptiveMeshN::FixArray(void) {
    //Set the fixed part of the array
    for (int i=0;i<num_dimensions;i++) {
	if (((pattern >> i) & 0x1) == 0x0) {
	    array[i]=x_min[i];
	}
    }

}



void PAdaptiveMeshN::Divide(Int_t num, Int_t my_layer) {
    if (is_divided) return; //already divided

    // cout << "Divide" << endl;

    is_divided=1;
    total_sub_size=1;
    layer=my_layer;


    if (num > MAX_SUBMESH_N) {
	Warning("Divide","num > MAX_SUBMESH_N");
	return;
    }

    //In a first step we try to figure out which of the
    //variable dimension have to be divided.
    //We loop over the dimensions and the max/min
    //of the "right" and "left" side of thy hypercube, respectively.

    Double_t max_atmax[MAX_DIMENSIONS],
	min_atmax[MAX_DIMENSIONS], // "left" side
	max_atmin[MAX_DIMENSIONS],
	min_atmin[MAX_DIMENSIONS]; // "right" side

    
    for (int i=0;i<num_dimensions;i++) {
	if (((pattern >> i) & 0x1) == 0x1) { 
	    //outer loop over variable dimensions
	    
	    //cout << "checking dim. " << i << endl;

	    //now try to get all possible min/max combinations of the remaining
	    //dimensions

	    Int_t other_variable_dimensions = variable_dimensions -1;

	    if (other_variable_dimensions) {
		Int_t num_combinations = (Int_t) pow(2.,other_variable_dimensions);
		for (int mm_pattern=0; mm_pattern<num_combinations; mm_pattern++) {
		    //next loop over the other variable dimensions
		    //for each dimension we take the min/max according
		    //to the mm_pattern
		    
		    Int_t pattern_position=0;
		    FixArray();
		    for (int j=0;j<num_dimensions;j++) {
			if ((((pattern >> j) & 0x1) == 0x1) && (i!=j)){
			    
			    //cout << "setting dim. " << j << endl;
			    Int_t local_bit=(mm_pattern>>pattern_position) & 0x1;
			    if (local_bit)
				array[j]=x_min[j];
			    else
				array[j]=x_max[j];
			    pattern_position++;
			}
		    }
		    //now evaluate the value
		    array[i]=x_max[i];
		    Double_t local_max=model->GetWeight(array);
		    array[i]=x_min[i];
		    Double_t local_min=model->GetWeight(array);
		    if (mm_pattern) {
			if (local_max > max_atmax[i])
			    max_atmax[i]=local_max;
			if (local_max < min_atmax[i])
			    min_atmax[i]=local_max;
			if (local_min > max_atmin[i])
			    max_atmin[i]=local_min;
			if (local_min < min_atmin[i])
			    min_atmin[i]=local_min;
		    } else {
			//first row in the table
			max_atmax[i]=min_atmax[i]=local_max;
			max_atmin[i]=min_atmin[i]=local_min;
		    }

		}

	    } else {
		//nothing left
		FixArray();
		array[i]=x_max[i];
		max_atmax[i]=min_atmax[i]=model->GetWeight(array);
		FixArray();
		array[i]=x_min[i];
		max_atmin[i]=min_atmin[i]=model->GetWeight(array);

	    }

	    if (max_atmax[i]<y_max) max_atmax[i]=y_max;

	    //Now we try to get an answer if the current dimension should be divided at all
	    Double_t diff1 = 0.;
	    if (min_atmin[i]) diff1 = max_atmax[i]/min_atmin[i];
	    Double_t abs_diff1 = fabs(max_atmax[i]-min_atmin[i]);
	    Double_t diff2 = 0.;
	    if (max_atmin[i]) diff2 = min_atmax[i]/max_atmin[i];
	    Double_t abs_diff2 = fabs(min_atmax[i]-max_atmin[i]);

 // 	    cout << diff1 << ":" << diff2 << ":" << abs_diff1
//  		 << ":" << abs_diff2 << endl;

	    if ((((diff1 > threshold_diff) || ((1./diff1)>threshold_diff)) &&
		 (abs_diff1>threshold_abs)) ||
		(((diff2 > threshold_diff) || ((1./diff2)>threshold_diff)) &&
		 (abs_diff2>threshold_abs))
		) {
		sub_size[i]=num;
		total_sub_size*=num;
	    }
	}//END outer loop over variable dimensions
    }

    //After this has been clarified the new sub-meshs have to created

    if (total_sub_size==1) { //forget it!
	total_sub_size=0;
	return;
    }

    sub_tree = new PAdaptiveMeshN[total_sub_size];

    for (int i=0;i<total_sub_size;i++) {
	sub_tree[i].SetDefaults(pattern, num_dimensions, 
				model, y_max);
	//cout << "sub_tree created" << endl;
    }

    sub_area = new Double_t[total_sub_size];

    
    

    // now we loop

    Int_t current_mesh_pos=0, current_dimension=0;

    Int_t sub_mesh_pos[MAX_DIMENSIONS];
    Double_t min[MAX_DIMENSIONS],max[MAX_DIMENSIONS];
    for (int i=0;i<num_dimensions;i++) sub_mesh_pos[i]=0;

    while (current_mesh_pos < total_sub_size) {
	//cout << current_mesh_pos<< ":"<< total_sub_size << endl;
	for (int i=0;i<num_dimensions;i++) {
	    //for each node we get the size of the sub-mesh
	    if (sub_size[i]>1) {
		min[i]=x_min[i]+((x_max[i]- x_min[i])/(Double_t)sub_size[i])*(Double_t)sub_mesh_pos[i];
		max[i]=min[i]+((x_max[i]- x_min[i])/(Double_t)sub_size[i]);
	    } else {
		min[i]=x_min[i];
		max[i]=x_max[i];
	    }

	    sub_tree[current_mesh_pos].SetRange(i,min[i],max[i]);
	    
	}

	//sub_tree[current_mesh_pos].PrintMesh();
	sub_tree[current_mesh_pos].SetThreshold(threshold_diff, threshold_abs);
	sub_tree[current_mesh_pos].SetMCPoints(mcpoints);
	sub_tree[current_mesh_pos].ReCalcYMax();
	sub_tree[current_mesh_pos].ReCalc();

	current_mesh_pos++;
	
	sub_mesh_pos[current_dimension]++;
	
	if (current_mesh_pos < total_sub_size) 
	    while(sub_mesh_pos[current_dimension] >= sub_size[current_dimension]) {
		
		for (int j=0;j<=current_dimension;j++) 
		    sub_mesh_pos[j]=0;
		
		current_dimension++;
		sub_mesh_pos[current_dimension]++;
	    }
	
	current_dimension=0;

    }

    ReCalcYMax();
    
    if (layer > 0) {
	//cout << "sub-divide" << endl;
	for (int i=0;i<total_sub_size;i++) {
	    sub_tree[i].Divide(num,layer-1);
	    sub_tree[i].ReCalcYMax();
	}
    }

    ReCalc();
    
}

void PAdaptiveMeshN::PrintMesh(void) {
    cout << "layer: " << layer <<  " area_size: " << area_size << endl;
    for (int j=0;j<num_dimensions;j++) {
	if ((((pattern >> j) & 0x1) == 0x1)) {
	    cout << " x_min[" << j << "]=" 
		 << x_min[j] 
		 << " x_max[" << j << "]=" 
		 << x_max[j];
	}
    }
    cout << " y_max=" << y_max << endl;

    for (int i=0;i<total_sub_size;i++) 
	sub_tree[i].PrintMesh();

}

void PAdaptiveMeshN::ReCalcYMax(void) {
    y_max=0.;
    Int_t num_combinations = (Int_t) pow(2.,variable_dimensions);

    for (int mm_pattern=0; mm_pattern<num_combinations; mm_pattern++) {
	//next loop over the other variable dimensions
	//for each dimension we take the min/max according
	//to the mm_pattern

	Int_t pattern_position=0;
	FixArray();
	for (int j=0;j<num_dimensions;j++) {
	    if ((((pattern >> j) & 0x1) == 0x1)){
		Int_t local_bit=(mm_pattern>>pattern_position) & 0x1;
		if (local_bit)
		    array[j]=x_min[j];
		else
		    array[j]=x_max[j];
		pattern_position++;
	    }
	}
	//now evaluate the value
	Double_t local_max=model->GetWeight(array);
	if (local_max>y_max) y_max = local_max;
    }

    for (int i=0; i<mcpoints; i++) {
	//In addition set some points in the middle
	FixArray();
	for (int j=0;j<num_dimensions;j++) {
	    if ((((pattern >> j) & 0x1) == 0x1)){
		array[j]=PUtils::sampleFlat()*(
		    GetXMax(j)-GetXMin(j)
		) + GetXMin(j);
	    }
	}
	Double_t local_max=model->GetWeight(array);
	if (local_max>y_max) y_max = local_max;
    }

    cout << "y_max" << y_max << endl;

}


// void PAdaptiveMesh::Draw(const Option_t*) {
//     if (sub_size) {
// 	for (int i=0;i<sub_size;i++) {
// 	    sub_tree[i]->Draw();
// 	}
	
//     } else {
// 	if (!line)
// 	    line = new TLine(x_min, y_max, x_max, y_max);
// 	else {
// 	    line->SetX1(x_min);
// 	    line->SetX2(x_max);
// 	    line->SetY1(x_max);
// 	    line->SetY2(x_max);
// 	}
// 	line->Draw();
//     }
//     tf1->Draw();
// }
