////////////////////////////////////////////////////////
//  Value container implementation file
//
//  Just a small container class for user-defined values  
//
////////////////////////////////////////////////////////


#include "PValues.h"
#include <iostream>
#include <string.h>

using namespace std;

PValues::PValues() {
    pointer=0;

}

PValues::PValues(const PValues  & p) {
    pointer=p.pointer;

    for (int i=0;i<pointer;i++) {
	array_id[i]  = p.array_id[i];
	array_val[i] = p.array_val[i];
    }
}

bool PValues::SetValue(int id , double val) {

    for (int i=0;i<pointer;i++) {
	if (array_id[i]==id) {
	    array_val[i]=val;
	    return kTRUE;
	}
    }
    if (pointer==MAX_VALUES) {
	return kFALSE;
    }

    array_val[pointer]=val;
    array_id[pointer]=id;

    pointer++;

    return kTRUE;

}

bool PValues::GetValue(int id , double * val) {
    for (int i=0;i<pointer;i++) {
	if (array_id[i]==id) {
	    *val=array_val[i];
	    return kTRUE;
	}
    }
    return kFALSE;

}

int  PValues::StringToValueID(char * st) {
    if (!strcmp(st,"t")) return T_MATRIX;
    if (!strcmp(st,"u")) return U_MATRIX;
    if (!strcmp(st,"tu")) return TU_MATRIX;
    return -1;
}

void PValues::Print(const Option_t* delme) const {
    for (int i=0;i<pointer;i++) {
	cout << "Value #"<< array_id[i] <<" is: " << array_val[i] << endl;
    }
}

 
ClassImp(PValues)
