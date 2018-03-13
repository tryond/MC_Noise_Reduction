#include "Sphere.h"

Sphere::Sphere(double rad_, Vec p_, Vec e_, const BRDF &brdf_) :
    rad(rad_), p(p_), e(e_), brdf(brdf_) {}

double Sphere::intersect(const Ray &r) const { // returns distance, 0 if nohit
    Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    double t, eps = 1e-4, b = op.dot(r.d), det = b*b-op.dot(op)+rad*rad;
    if ( det<0 ) return 0; else det = sqrt(det);
    return (t = b-det)>eps ? t : ((t = b+det)>eps ? t : 0);
}
