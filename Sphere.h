#ifndef SPHERE_H
#define SPHERE_H

#include <cmath>
#include "Vec.h"
#include "BRDF.h"

struct Sphere {
    Vec p, e;           // position, emitted radiance
    double rad;         // radius
    const BRDF &brdf;   // BRDF

    Sphere(double rad_, Vec p_, Vec e_, const BRDF &brdf_);

    double intersect(const Ray &r) const;
};

bool intersect(const Ray &r, double &t, int &id) {
    double n = sizeof(spheres)/sizeof(Sphere), d, inf = t = 1e20;
    for ( int i = int(n); i--;) if ( (d = spheres[i].intersect(r))&&d<t ) { t = d; id = i; }
    return t<inf;
}

#endif