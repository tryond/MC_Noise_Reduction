// Submitter: tryond(tryon,daniel) - 20621204

// Final Project

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <random>
#include <vector>
#include <omp.h>
#include <iostream>

// *******************************************************************
#include <time.h>
#include <algorithm>
#include <ctime>

#include "RNG.h"
#include "Vec.h"

#include "SimpleFilter.h"
#include "ObjectFilter.h"
#include "BrdfFilter.h"
#include "NormalFilter.h"
#include "DofFilter.h"
// *******************************************************************

#define PI 3.1415926535897932384626433832795

// Received Radiance Values
const int rrDepth = 5;
const double survivalProbability = 0.9;

// Thread-safe random number generator
RNG rng;

struct BRDF {
    virtual Vec eval(const Vec &n, const Vec &o, const Vec &i) const = 0;
    virtual void sample(const Vec &n, const Vec &o, Vec &i, double &pdf) const = 0;
    virtual bool isSpecular() const = 0;
};


/*
 * Shapes
 */
struct Sphere {
    Vec p, e;           // position, emitted radiance
    double rad;         // radius
    const BRDF &brdf;   // BRDF

    Sphere(double rad_, Vec p_, Vec e_, const BRDF &brdf_) :
        rad(rad_), p(p_), e(e_), brdf(brdf_) {}

    double intersect(const Ray &r) const { // returns distance, 0 if nohit
        Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t, eps = 1e-4, b = op.dot(r.d), det = b*b-op.dot(op)+rad*rad;
        if ( det<0 ) return 0; else det = sqrt(det);
        return (t = b-det)>eps ? t : ((t = b+det)>eps ? t : 0);
    }
};



/*
 * Sampling functions
 */

inline void createLocalCoord(const Vec &n, Vec &u, Vec &v, Vec &w) {
    w = n;
    u = ((std::abs(w.x)>.1 ? Vec(0, 1) : Vec(1)).cross(w)).normalize();
    v = w.cross(u);
}

Vec uniformRandomPSA(const Vec &n) {
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


/*
 * BRDFs
 */
// Ideal diffuse BRDF
struct DiffuseBRDF : public BRDF {
    DiffuseBRDF(Vec kd_) : kd(kd_) {}

    Vec eval(const Vec &n, const Vec &o, const Vec &i) const {
        return kd * (1.0/PI);
    }

    void sample(const Vec &n, const Vec &o, Vec &i, double &pdf) const {

        // i is random unit vector in hemisphere around normal, n
        i = uniformRandomPSA(n);

        // pdf = cos(theta) / PI
        pdf = i.dot(n)/PI; // dot of two unit vectors is cos(theta)
    }

    bool isSpecular() const {
        return false;
    }

    Vec kd;
};

struct SpecularBRDF : public BRDF {
    SpecularBRDF(Vec ks_) : ks(ks_) {}
    
    Vec reflectVec(const Vec &n, const Vec &o) const {
        return (n*n.dot(o)*2.0) - o; 
    }

    Vec eval(const Vec &n, const Vec &o, const Vec &i) const {            
        Vec ref = reflectVec(n,o);
        if (ref.x != i.x || ref.y != i.y || ref.z != i.z) 
            return Vec();
        return ks * (1.0 / n.dot(i));
    }

    void sample(const Vec &n, const Vec &o, Vec &i, double &pdf) const {
        // Reflection Vector
        i = reflectVec(n,o);   // i = 2<o,n>n - o
        pdf = 1.0;
    }

    bool isSpecular() const {
        return true;
    }

    Vec ks;    
};


/*
 * Scene configuration
 */
// Pre-defined BRDFs
const DiffuseBRDF leftWall(Vec(.75,.25,.25)),
                  rightWall(Vec(.25,.25,.75)),
                  otherWall(Vec(.75,.75,.75)),
                  blackSurf(Vec(0.0,0.0,0.0)),
                  brightSurf(Vec(0.9,0.9,0.9));

const SpecularBRDF mirrorSurf(Vec(0.999,0.999,0.999));

// Scene: list of spheres
const Sphere spheres[] = {
    Sphere(1e5,  Vec(1e5+1,40.8,81.6),   Vec(),         leftWall),   // Left
    Sphere(1e5,  Vec(-1e5+99,40.8,81.6), Vec(),         rightWall),  // Right
    Sphere(1e5,  Vec(50,40.8, 1e5),      Vec(),         otherWall),  // Back
    Sphere(1e5,  Vec(50, 1e5, 81.6),     Vec(),         otherWall),  // Bottom
    Sphere(1e5,  Vec(50,-1e5+81.6,81.6), Vec(),         otherWall),  // Top
    Sphere(16.5, Vec(27,16.5,47),        Vec(),         brightSurf), // Ball 1
    Sphere(16.5, Vec(73,16.5,78),        Vec(),         mirrorSurf), // Ball 2
    Sphere(5.0,  Vec(50,70.0,81.6),      Vec(50,50,50), blackSurf)   // Light
};


/*
 * Global functions
 */
bool intersect(const Ray &r, double &t, int &id) {
    double n = sizeof(spheres)/sizeof(Sphere), d, inf = t = 1e20;
    for ( int i = int(n); i--;) if ( (d = spheres[i].intersect(r))&&d<t ) { t = d; id = i; }
    return t<inf;
}

Vec directRadiance(const Ray &ray, const Vec &xn, int idx, const BRDF &brdf) {

    Vec dirRad(0,0,0);                          // Direct = zero vec if not visible

    if (brdf.isSpecular()) {
        return dirRad;
    }    


    Vec x = ray.o;                              // Origin
    Vec o = ray.d;                              // Direction
    
    Vec y, yn;
    double ypdf;
    const Sphere &lum = spheres[7];
    luminaireSample(lum, y, yn, ypdf);          // Get pos on lum, normal, and pdf
    Vec xin = (y-x).normalize();                // Direct Light on x from y
    double r2 = (y-x).dot(y-x);                 // Magnitude of y-x squared
    
    Vec yin = xin*-1.0;                         // Vector from x to y 
    if (yn.dot(yin) < 0) 
        yn = yn*-1.0;                           // Normal correction


    bool vis;                                   // True if x visible to y
    int idx2 = 0;                               // id at point x
    int idy2 = 0; 
    double t;                                   // id at point y
    if (!intersect(Ray(y,yin),t,idx2) || 
        !intersect(Ray(x,xin),t,idy2)) {
        vis = false;
    }
    else if (idx2 == idx && idy2 == 7) 
        vis = true;
    else 
        vis = false;


    if (vis) {
        double xdot = xn.dot(xin);
        double ydot = yn.dot(yin);
        double rdot = r2*ypdf;
        dirRad = lum.e.mult(brdf.eval(xn, o, xin))*(xdot*ydot/rdot);
    }

    return dirRad;
}

Vec reflectedRadiance(const Ray &r, int depth){
    double t;                                   // Distance to intersection
    
    // *******************************************************************
    // should this be -1?
    int id = 0;                                 // id of intersected sphere
    
    if (!intersect(r, t, id)) 
        return Vec();   // if miss, return black
    const Sphere &obj = spheres[id];            // the hit object
    
    Vec x = r.o + r.d*t;                        // The intersection point
    Vec o = (Vec() - r.d).normalize();          // The outgoing direction (= -r.d)
    
    Vec xn = (x - obj.p).normalize();            // The normal direction
    if (xn.dot(o) < 0) 
        xn = xn*-1.0;

    const BRDF &brdf = obj.brdf;                // Surface BRDF at x

    Vec retRad = directRadiance(Ray(x,o), xn, id, brdf);

    // INDIRECT
    double p = 1.0;
    if (depth > rrDepth)
        p = survivalProbability;
    
    if (rng() < p){
        Vec xin2;
        double xpdf;
        brdf.sample(xn, o, xin2, xpdf);
        Ray rx(x, xin2);
        retRad = retRad + brdf.eval(xn, o, xin2).mult(reflectedRadiance(rx, depth+1))*(xn.dot(xin2)/(xpdf*p));
    }
    return retRad;
}

Vec receivedRadiance(const Ray &r, int depth, bool flag) {
    
    double t; 
    int id = 0;                                 	// id of intersected sphere
    
    if (!intersect(r, t, id)) 
        return Vec();                           // if miss, return black
    const Sphere &obj = spheres[id];            // the hit object
  
    Vec x = r.o + r.d*t;                        // The intersection point
    Vec o = (Vec() - r.d).normalize();          // The outgoing direction (= -r.d)

    Vec n = (x - obj.p).normalize();            // The normal direction
    if ( n.dot(o) < 0 ) n = n*-1.0;

    Vec rad = obj.e;                            // Emitted radiance
    
    if (!obj.brdf.isSpecular())
        return rad + reflectedRadiance(r, depth);   // Emitted plus reflected vector

    else {

        double p = 1.0;
        if (depth > rrDepth)
            p = survivalProbability;

        if (rng() < p) {

            const BRDF &brdf = obj.brdf;
            double pdf;                             // PDF for brdf
            Vec i;                                  // BRDF sample incoming

            brdf.sample(n,o,i,pdf);                 // Get incoming and pdf

            // Find y, and -w
            Ray r2(x,i);

            // Le + Lr = Le + Integrate(Incoming * BRDF * cos(theta))
            return rad + receivedRadiance(r2,depth+1,flag).mult(brdf.eval(n,o,i)) * (n.dot(i)/(pdf*p));
        }
    }
    return rad;


}

// Camera position & direction
// const Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).normalize());
const Ray cam(Vec(50, 52, 200.0), Vec(0, -0.042612, -1).normalize());


/*
 * Main function (do not modify)
 */
int main(int argc, char *argv[]) {
    
    // setup elapsed time
    struct timespec start, finish;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &start);

    int nworkers = omp_get_num_procs();
    omp_set_num_threads(nworkers);
    rng.init(nworkers);

    // int w = 480, h = 360, samps = argc==2 ? atoi(argv[1])/4 : 1; // # samples
    int w = 480, h = 360, samps = argc >= 2 ? atoi(argv[1])/4 : 1; // # samples

    // *******************************************************************
    int sMean = argc >= 3 ? atoi(argv[2]) : 0;
    int sMed = argc >= 4 ? atoi(argv[3]) : 0;
    int sGaus = argc >= 5 ? atoi(argv[4]) : 0;

    int oMean = argc >= 6 ? atoi(argv[5]) : 0;
    int oMed = argc >= 7 ? atoi(argv[6]) : 0;
    int oGaus = argc >= 8 ? atoi(argv[7]) : 0;

    int bMean = argc >= 9 ? atoi(argv[8]) : 0;
    int bGaus = argc >= 10 ? atoi(argv[9]) : 0;

    int nMean = argc >= 11 ? atoi(argv[10]) : 0;
    int nGaus = argc >= 12 ? atoi(argv[11]) : 0;

    int dMean = argc >= 13 ? atoi(argv[12]) : 0;
    int dFoc = argc >= 14 ? atoi(argv[13]) : 0;
    // *******************************************************************

    Vec cx = Vec(w*.5135/h), cy = (cx.cross(cam.d)).normalize()*.5135;
    std::vector<Vec> c(w*h);
    
    // *******************************************************************
    
    // pixel's object
    std::vector<int> o(w*h);

    // pixel's reflection vector (-1 if not specular)
    std::vector<int> r(w*h);

    // pixel's normal vector
    std::vector<Vec> n(w*h);

    // pixel's depth
    std::vector<int> d(w*h);

    // object id 
    int id;

#pragma omp parallel for schedule(dynamic, 1)
    for ( int y = 0; y < h; y++ ) {
        for ( int x = 0; x < w; x++ ) {
           
            const int i = (h - y - 1)*w + x;

            // *******************************************************************
            // COMPUTE OBJECT HERE
            double dist;
            id = 0;
            // Shoot ray through the center of pixel i
            Vec dir = cx*(((x+.5)/w) - .5) + cy*(((y+.5)/h) - .5) + cam.d;
            Ray r1(cam.o,dir.normalize());
            intersect(r1,dist,id);
            o[i] = id;

            // *******************************************************************
            // ASSIGN PIXEL DEPTH
            d[i] = dist;

            // *******************************************************************
            // COMPUTE PIXEL NORMAL
            Vec orig = r1.o + r1.d*dist;                    // intersection point
            Vec out = (Vec()-r1.d).normalize();              // The outgoing direction (= -r.d)
            Vec norm = (orig-spheres[id].p).normalize();      // The normal direction 
            if (norm.dot(out) < 0) 
                norm = norm*-1.0;
            n[i] = norm;

            // *******************************************************************
            // COMPUTE REFLECT OBJECT HERE
            if (spheres[id].brdf.isSpecular()) {
                Vec refVec = (norm*norm.dot(orig)*2.0) - orig;
                Ray refRay(orig,refVec.normalize());
                if (intersect(refRay,dist,id))
                    r[i] = id;
                else
                    r[i] = -1;
            }
            else {
                r[i] = -1;
            }

            for ( int sy = 0; sy < 2; ++sy ) {
                for ( int sx = 0; sx < 2; ++sx ) {
                    Vec rad;
                    int id_sum = 0;
                    for ( int s = 0; s<samps; s++ ) {
                        double r1 = 2*rng(), dx = r1<1 ? sqrt(r1)-1 : 1-sqrt(2-r1);
                        double r2 = 2*rng(), dy = r2<1 ? sqrt(r2)-1 : 1-sqrt(2-r2);
                        Vec d = cx*(((sx+.5 + dx)/2 + x)/w - .5) + cy*(((sy+.5 + dy)/2 + y)/h - .5) + cam.d;
                        rad = rad + receivedRadiance(Ray(cam.o, d.normalize()), 1, true)*(1./samps);
                    }
                    c[i] = c[i] + Vec(clamp(rad.x), clamp(rad.y), clamp(rad.z))*.25;
                }
            }
            
        }
#pragma omp critical
        fprintf(stderr,"\rRendering (%d spp) %6.2f%%",samps*4,100.*y/(h-1));
    }
    fprintf(stderr, "\n");

    // *******************************************************************
    // APPLY FILTERS

    // Simple Filters
    c = simpleMean(c,h,w,sMean);            // arg 1
    c = simpleMedian(c,h,w,sMed);           // arg 2
    c = simpleGaussian(c,h,w,sGaus);        // arg 3

    // Object Filters
    c = objectMean(c,o,h,w,oMean);          // arg 4
    c = objectMedian(c,o,h,w,oMed);         // arg 5
    c = objectGaussian(c,o,h,w,oGaus);      // arg 6

    // BRDF Filters
    c = brdfMean(c,o,r,h,w,bMean);          // arg 7
    c = brdfGaussian(c,o,r,h,w,bGaus);      // arg 8

    // Normal Filters
    c = normalMean(c,o,n,h,w,nMean);        // arg 9
    c = normalGaussian(c,o,n,h,w,nGaus);    // arg 10

    // DOF Filter
    c = dofMean(c,o,d,h,w,dMean,dFoc);      // arg 11, arg 12
    // *******************************************************************


    // Write resulting image to a PPM file
    FILE *f = fopen("image.ppm", "w");
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for ( int i = 0; i<w*h; i++ )
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
    fclose(f);
	
	// Get time elapsed
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    std::cout << "Time: " << elapsed << "s\n";

    return 0;
}
