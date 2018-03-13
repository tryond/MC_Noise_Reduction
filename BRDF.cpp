#include "BRDF.h"


DiffuseBRDF::DiffuseBRDF(Vec kd_) : kd(kd_) {}

Vec DiffuseBRDF::eval(const Vec &n, const Vec &o, const Vec &i) const {
    return kd * (1.0/PI);
}

void DiffuseBRDF::sample(const Vec &n, const Vec &o, Vec &i, double &pdf) const {
    // i is random unit vector in hemisphere around normal, n
    i = uniformRandomPSA(n);
    pdf = i.dot(n)/PI; // dot of two unit vectors is cos(theta)
}

bool DiffuseBRDF::isSpecular() const {
    return false;
}



    
SpecularBRDF::SpecularBRDF(Vec ks_) : ks(ks_) {}

Vec SpecularBRDF::reflectVec(const Vec &n, const Vec &o) const {
    return (n*n.dot(o)*2.0) - o; 
}

Vec SpecularBRDF::eval(const Vec &n, const Vec &o, const Vec &i) const {            
    Vec ref = reflectVec(n,o);
    if (ref.x != i.x || ref.y != i.y || ref.z != i.z) 
        return Vec();
    return ks * (1.0 / n.dot(i));
}

void SpecularBRDF::sample(const Vec &n, const Vec &o, Vec &i, double &pdf) const {
    // Reflection Vector
    i = reflectVec(n,o);   // i = 2<o,n>n - o
    pdf = 1.0;
}

bool SpecularBRDF::isSpecular() const {
    return true;
}

