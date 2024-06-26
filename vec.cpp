#include "vec.h"

vec scaleproduct (double s, const vec &A){
    vec res;
    res.x = s*A.x;
    res.y = s*A.y;
    res.z = s*A.z;
    return res;
}
double dotproduct(const vec &A,const vec &B){
    return A.x*B.x + A.y*B.y + A.z*B.z;
}

double veclength(const vec &A){
    return sqrt(A.x*A.x + A.y*A.y + A.z*A.z);
}

vec crossproduct(const vec &A, const vec &B){
    vec res;
    res.x = A.y*B.z - A.z*B.y;
    res.y = A.z*B.x - A.x*B.z;
    res.z = A.x*B.y - A.y*B.x;
    return res;
}

vec diff_vec(const vec &A, const vec &B){
    vec res;
    res.x = A.x-B.x;
    res.y = A.y-B.y;
    res.z = A.z-B.z;
    return res;
}

double projlenJonR(const vec &J,const vec &R){
    return veclength(J)* dotproduct(J,R)/(veclength(J)*veclength(R));
}

vec projvecJonR(const vec &J, const vec &R){
    return scaleproduct((dotproduct(J,R))/(veclength(R)*veclength(R)), R);
}

vec average_vectors(vector<vec> &vecs){
    vec res;
    double x,y,z;
    x=y=z=0;
    for(int i=0; i< vecs.size(); ++i){
        x += vecs[i].x;
        y += vecs[i].y;
        z += vecs[i].z;
    }
    res.x = x/vecs.size();
    res.y = y/vecs.size();
    res.z = z/vecs.size();
    return res;
}


vec sum_vec(const vec &A, const vec &B){
    vec res;
    res.x = A.x+B.x;
    res.y = A.y+B.y;
    res.z = A.z+B.z;
    return res;
}

vec unit_vec(const vec &A){
    vec res;
    res.x = A.x/veclength(A);
    res.y = A.y/veclength(A);
    res.z = A.z/veclength(A);
    return res;
}