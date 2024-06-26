#ifndef vec_H_INCLUDE
#define vec_H_INCLUDE
#include<math.h>
#include<vector>
#include <sstream>
#include <iomanip>
using namespace std;
struct vec
{
    double x,y,z;
    void print(ostream& stream) const{
		stream<<setprecision(11) << fixed << right << setw(20) <<this->x
			  <<setprecision(11) << fixed << right << setw(20) <<this->y 
			  <<setprecision(11) << fixed << right << setw(20) <<this->z <<endl;
	}
};
vec diff_vec(const vec &A, const vec &B);
vec unit_vec(const vec &A);
vec sum_vec(const vec &A, const vec &B);
double projlenJonR(const vec &J, const vec &R);
vec scaleproduct (double s, const vec &A);
double dotproduct(const vec &A, const vec &B);
double veclength(const vec &A);
vec crossproduct(const vec &A, const vec &B);
vec projvecJonR(const vec &J, const vec &R);
vec average_vectors(vector<vec> &vecs);

#endif