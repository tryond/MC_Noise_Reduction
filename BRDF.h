#ifndef BRDF_H
#define BRDF_H

#include "Vec.h"
#include "RNG.h"
#include "Sphere.h"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

const double PI = 3.1415926535897932384626433832795;

// RNG rng;

struct BRDF {
    virtual Vec eval(const Vec &n, const Vec &o, const Vec &i) const = 0;
    virtual void sample(const Vec &n, const Vec &o, Vec &i, double &pdf) const = 0;
    virtual bool isSpecular() const = 0;
};

// Ideal diffuse BRDF
struct DiffuseBRDF : public BRDF {
    
    DiffuseBRDF(Vec kd_);

    Vec eval(const Vec &n, const Vec &o, const Vec &i) const;
    
    void sample(const Vec &n, const Vec &o, Vec &i, double &pdf) const;

    bool isSpecular() const;

    Vec kd;
};

struct SpecularBRDF : public BRDF {
    
    SpecularBRDF(Vec ks_);
    
    Vec reflectVec(const Vec &n, const Vec &o) const;

    Vec eval(const Vec &n, const Vec &o, const Vec &i) const;

    void sample(const Vec &n, const Vec &o, Vec &i, double &pdf) const;

    bool isSpecular() const;

    Vec ks;    
};

inline void createLocalCoord(const Vec &n, Vec &u, Vec &v, Vec &w) {
    w = n;
    u = ((std::abs(w.x)>.1 ? Vec(0, 1) : Vec(1)).cross(w)).normalize();
    v = w.cross(u);
}

Vec uniformRandomPSA(const Vec &n) {
    RNG rng;

    double z = sqrt(rng());
    double r = sqrt(1.0 - pow(z,2));
    double phi = 2.0 * PI * rng();
    double x = r * cos(phi);
    double y = r * sin(phi);

    Vec u, v, w;
    createLocalCoord(n, u, v, w);

    // return unit vector in hemisphere around n
    return u*x + v*y + w*z;
}

void luminaireSample(const Sphere &lum, Vec &y, Vec &n, double &pdf){
    
    RNG rng;

    const Vec c = lum.p;
    const double r = lum.rad;

    double zp = 2.0*rng() - 1.0;

    double rand = rng();
    double xp = sqrt(1.0-pow(zp,2))*cos(2.0*PI*rand);
    double yp = sqrt(1.0-pow(zp,2))*sin(2.0*PI*rand);
    
    pdf = 1.0/(4.0*PI*pow(r,2));  // 1.0/Surface Area
    n = Vec(xp,yp,zp).normalize();   // normal direction

    // y is point on surface of luminaire
    y = c + n*r;  // y = (x0,y0,z0) + r0(xp,yp,zp)
}

#endif