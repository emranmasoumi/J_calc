#ifndef read_geom_H_INCLUDE
#define read_geom_H_INCLUDE
#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include "vec.h"

using namespace std;

struct FluxJ{
    string IELEMENT;
    string IIELEMENT;
    vec J;
    void print(ostream& stream) const{
		stream<<this->IELEMENT<<"\t"<<this->IIELEMENT<<"\t";
		stream<<setprecision(11) << fixed << right << setw(20) <<this->J.x
			  <<setprecision(11) << fixed << right << setw(20) <<this->J.y 
			  <<setprecision(11) << fixed << right << setw(20) <<this->J.z <<endl;
	}
};

struct Atom{
    vec pos;
    double charge;
    string ELEMENT;

    void print(ostream& stream) const{
		stream<<this->ELEMENT<<"\t";
		stream<<setprecision(11) << fixed << right << setw(15) <<this->pos.x
			  <<setprecision(11) << fixed << right << setw(15) <<this->pos.y 
			  <<setprecision(11) << fixed << right << setw(15) <<this->pos.z <<endl;
	}
};



void read_geometry(vector<Atom> &Atoms, string filename);

void read_extracted_FluxJ(vector<FluxJ> &FluxJs, string filename);

#endif