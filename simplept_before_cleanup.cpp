// Submitter: tryond(tryon,daniel) - 20621204

// Task 3-2

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <random>
#include <vector>
#include <omp.h>
#include <iostream>

#include <time.h>
// *******************************************************************
#include <algorithm>
#include <ctime>

// *******************************************************************
#include "RNG.h"
#include "Vec.h"



#define PI 3.1415926535897932384626433832795

// Received Radiance Values
const int rrDepth = 5;
const double survivalProbability = 0.9;


// Thread-safe random number generator
RNG rng;


/*
 * Thread-safe random number generator
 */

/*
struct RNG {
    RNG() : distrb(0.0, 1.0), engines() {}

    void init(int nworkers) {
        std::random_device rd;
        engines.resize(nworkers);
        for ( int i = 0; i < nworkers; ++i )
            engines[i].seed(rd());
    }

    double operator()() {
        int id = omp_get_thread_num();
        return distrb(engines[id]);
    }

    std::uniform_real_distribution<double> distrb;
    std::vector<std::mt19937> engines;
} rng;
*/



/*
 * Basic data types
 */

/*
const double epsilon = 0.0001;

struct Vec {
    double x, y, z;


// ****************************************************************
	double mag; // be sure to check initializer
	bool operator< (const Vec &b) const { return mag < b.mag; } 
	bool operator== (const Vec &b) const { return abs(x-b.x) > epsilon ? 0 : 
												  abs(y-b.y) > epsilon ? 0 : 
												  abs(z-b.z) > epsilon ? 0 : 1; }
	Vec& operator= (const Vec &b) {
		x = b.x;
		y = b.y;
		z = b.z;
		mag = b.mag;
		return *this;
	}
	
	
	Vec average(const Vec &b) const { return Vec((x+b.x)/2.0, (y+b.y)/2.0, (z+b.z)/2.0); }	
	
	Vec(double x_ = 0, double y_ = 0, double z_ = 0) { x = x_; y = y_; z = z_; mag = x_*x_ + y_*y_ + z_*z_; }
	
	Vec operator/ (double b) const { return Vec(x/b, y/b, z/b); } 	
	
	
// ****************************************************************

    Vec operator+ (const Vec &b) const  { return Vec(x+b.x, y+b.y, z+b.z); }
    Vec operator- (const Vec &b) const  { return Vec(x-b.x, y-b.y, z-b.z); }
    Vec operator* (double b) const      { return Vec(x*b, y*b, z*b); }

    Vec mult(const Vec &b) const        { return Vec(x*b.x, y*b.y, z*b.z); }
    Vec& normalize()                    { return *this = *this * (1.0/std::sqrt(x*x+y*y+z*z)); }
    double dot(const Vec &b) const      { return x*b.x+y*b.y+z*b.z; }
    Vec cross(const Vec &b) const        { return Vec(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x); }
  
};

std::ostream& operator<<(std::ostream& out, const Vec &b)
{
    out << "(" << b.x << ", " << b.y << ", " << b.z << ")" << " = " << b.mag;
    return out;
}
*/

struct Ray {
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

struct BRDF {
    virtual Vec eval(const Vec &n, const Vec &o, const Vec &i) const = 0;
    virtual void sample(const Vec &n, const Vec &o, Vec &i, double &pdf) const = 0;
    virtual bool isSpecular() const = 0;
};


/*
 * Utility functions
 */

inline double clamp(double x)   {
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

inline int toInt(double x) {
    return static_cast<int>(std::pow(clamp(x), 1.0/2.2)*255+.5);
}


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

// Camera position & direction
const Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).normalize());


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

// ****************************************************
// added "int &id" to definition

Vec receivedRadiance(const Ray &r, int depth, bool flag) {
    
    double t; 
    // ****************************************************                                
    // Distance to intersection
    // int id = 0;
    // should this be -1?
    int id = 0;                                 	// id of intersected sphere
    
    if (!intersect(r, t, id)) 
        return Vec();                           // if miss, return black
    const Sphere &obj = spheres[id];            // the hit object
  

    // HERE
    Vec x = r.o + r.d*t;                        // The intersection point
    Vec o = (Vec() - r.d).normalize();          // The outgoing direction (= -r.d)

    Vec n = (x - obj.p).normalize();            // The normal direction
    if ( n.dot(o) < 0 ) n = n*-1.0;
    // HERE




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
            
            // *******************************************************************
            // int id2 = id;
            return rad + receivedRadiance(r2,depth+1,flag).mult(brdf.eval(n,o,i)) * (n.dot(i)/(pdf*p));
        }
    }
    return rad;


}

// *******************************************************************
// *******************************************************************
// *******************************************************************
std::vector<Vec> boxObjFilter(std::vector<Vec> c, std::vector<int> o, int h, int w, int size) {

    if (size <= 1)
        return c;

    // size of filter should be odd
    if (size % 2 == 0) {
        ++size;
    }

    // total size of image
    int pix = h * w;
    
    // initialize return vector
    std::vector<Vec> b(pix);
    
    // corner move
    int move = (size-1)/2;
    
    // corner index
    int ci, cx, cy;
    
    // point indices
    int x, y;
        
    // potential x, y, and i
    int px, py, pi;
        
    // filter pixel sum
    Vec pixelSum;
    int pixelNum;
        
    // For each pixel in the source
    for (int i = 0; i < pix; ++i) {
        
        // Set Vec to center pixel
        pixelSum = c[i];
        pixelNum = 1;
        
        // Find x and y values
        x = i % w;
        y = ( (h*w) - i - 1) / w;
        
        // Find corner indices
        cx = x - move;
        cy = y - move;
    
        for (int v = 0; v < size; ++v) {
            for (int u = 0; u < size; ++u) {
            
                // Get potential pixel's x and y coords
                px = u + cx;
                py = v + cy;
                
                // DEBUG
                if (i == 0) {
                        // fprintf(stderr, "potential (%i,%i): ", px, py);
                }


                // Check that index is valid
                if (px >= 0 && px < w && py >= 0 && py < h) {
                    
                    // Compute potential pixel's index
                    pi = (h - py - 1) * w + px;
                    
                    // DEBUG
                    pixelSum = pixelSum + c[pi];
                    ++pixelNum;

                    
                    // Check that object id's are the same
                    /*
                    if (o[pi] == o[i]) {
                        pixelSum = pixelSum + c[pi];
                        ++pixelNum;
                    }
                    */
                    
                }           
            }
        } // filterPix has all contributing pixel values
    
        // Compute the median value
        b[i] = pixelSum * (1.0/(double)pixelNum);
        // b[i] = pixelSum;
    }
    
    // Return the vector of median values
    return b;           
}









Vec medianVec(std::vector<Vec> v) {

    sort(v.begin(), v.end());
    
    if (v.size() % 2 == 0)
        return v[v.size()/2].average(v[v.size()/2-1]);

    return v[v.size()/2];
}

std::vector<Vec> medObjFilter(std::vector<Vec> c, std::vector<int> o, int h, int w, int size) {

    if (size <= 1)
        return c;

    // size of filter should be odd
    if (size % 2 == 0) {
        ++size;
    }

    // total size of image
    int pix = h * w;
    
    // initialize return vector
    std::vector<Vec> m(pix);
    
    // corner move
    int move = (size-1)/2;
    
    // corner index
    int ci, cx, cy;
    
    // point indices
    int x, y;
        
    // potential x, y, and i
    int px, py, pi;
        
    // For each pixel in the source
    for (int i = 0; i < pix; ++i) {
        
        // Set Vec to center pixel
        std::vector<Vec> filterPix;
        filterPix.push_back(c[i]);
        
        // Find x and y values
        x = i % w;
        y = ( (h*w) - i - 1) / w;
        
        // Find corner indices
        cx = x - move;
        cy = y - move;
    
        for (int v = 0; v < size; ++v) {
            for (int u = 0; u < size; ++u) {
            
                // Get potential pixel's x and y coords
                px = u + cx;
                py = v + cy;
                
                // Check that index is valid
                if (px >= 0 && px < w && py >= 0 && py < h) {
                    
                    // Compute potential pixel's index
                    pi = (h - py - 1) * w + px;
                    
                    // DEBUG
                    // filterPix.push_back(c[pi]);

                    
                    // Check that object id's are the same
                    if (o[pi] == o[i]) {
                        filterPix.push_back(c[pi]);
                    }
                    
                    
                }           
            }
        } // filterPix has all contributing pixel values
    
        // Compute the median value
        m[i] = medianVec(filterPix);
        // std::cout << o[i] << ' ';
        // Vec v((double)o[i]*(255.0/20.0),(double)o[i]*(255.0/20.0),(double)o[i]*(255.0/20.0));
        // m[i] = v;
    }
    
    // Return the vector of median values
    return m;           
}


// **************************************************************************************

int medianInt(std::vector<int> v) {

    sort(v.begin(), v.end());
    
    /*
    if (v.size() % 2 == 0)
        return v[v.size()/2] + v[v.size()/2-1];
    */

    return v[v.size()/2];
}


std::vector<int> medArrFilter(std::vector<int> o, int h, int w, int size) {

    if (size <= 1)
        return o;

    // size of filter should be odd
    if (size % 2 == 0) {
        ++size;
    }

    // total size of image
    int pix = h * w;
    
    // initialize return vector
    std::vector<int> m(pix);
    
    // corner move
    int move = (size-1)/2;
    
    // corner index
    int ci, cx, cy;
    
    // point indices
    int x, y;
        
    // potential x, y, and i
    int px, py, pi;
        
    // For each pixel in the source
    for (int i = 0; i < pix; ++i) {
        
        // Set Vec to center pixel
        std::vector<int> filterPix;
        filterPix.push_back(o[i]);
        
        // Find x and y values
        x = i % w;
        y = ( (h*w) - i - 1) / w;
        
        // Find corner indices
        cx = x - move;
        cy = y - move;
    
        for (int v = 0; v < size; ++v) {
            for (int u = 0; u < size; ++u) {
            
                // Get potential pixel's x and y coords
                px = u + cx;
                py = v + cy;
                
                // Check that index is valid
                if (px >= 0 && px < w && py >= 0 && py < h) {
                    
                    // Compute potential pixel's index
                    pi = (h - py - 1) * w + px;
                    
                    // DEBUG
                    filterPix.push_back(o[pi]);

                    // Check that object id's are the same
                    /*
                    if (o[pi] == o[i]) {
                        filterPix.push_back(c[pi]);
                    }
                    */
                }           
            }
        } // filterPix has all contributing pixel values
    
        // Compute the median value
        m[i] = medianInt(filterPix);
        // std::cout << o[i] << ' ';
        // Vec v((double)o[i]*(255.0/20.0),(double)o[i]*(255.0/20.0),(double)o[i]*(255.0/20.0));
        // m[i] = v;
    }
    
    // Return the vector of median values
    return m;           
}









std::vector<std::vector<int>> createGaussian(int size) {

    // Mid point (max at v[max][max])
    int max = size/2;
  
    //.initialize kernel
    std::vector<std::vector<int>> v(size,std::vector<int>(size));
    for (int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            v[i][j] = 0;
    
    // set first row (Pascal's Triangle)
    int x = 1;
    for (int i = 0; i <= max; ++i) {
        v[0][i] = x;
        x = x * ((size-1) - i) / (i + 1);
    }
    
    // fill gaussian  
    for (int j = 1; j <= max; ++j) 
        for (int i = j; i <= max; ++i)
            v[j][i] = v[0][i] * v[0][j];
            
    // flip diagonally in first quadrant
    for (int j = 0; j <= max; ++j) 
        for (int i = j; i <= max; ++i)             
            v[i][j] = v[j][i];
    
    // flip accross vertical line in upper half
    for (int j = 0; j <= max; ++j) 
        for (int i = 0; i < max; ++i) 
            v[j][size-i-1] = v[j][i];
    
    // flip accross horizontal line in lower half
    for (int i = 0; i < size; ++i) 
        for (int j = 0; j < max; ++j) 
            v[size-j-1][i] = v[j][i];
    
    // DEBUG
    /*
    for (int j = 0; j < size; ++j) {
        for (int i = 0; i < size; ++i) {
            std::cout << v[j][i] << ' ';       
        }        
        std::cout << std::endl;
    }
    */

    return v;
}

std::vector<Vec> gausObjFilter(std::vector<Vec> c, std::vector<int> o, int h, int w, int size) {

    if (size <= 1)
        return c;

    // size of filter should be odd
    if (size % 2 == 0) {
        ++size;
    }

    // create the gaussian filter
    std::vector<std::vector<int>> kernel = createGaussian(size);

    // total size of image
    int pix = h * w;
    
    // initialize return vector
    std::vector<Vec> g(pix);
    
    // corner move
    int move = (size-1)/2;
    
    // corner index
    int ci, cx, cy;
    
    // point indices
    int x, y;
        
    // potential x, y, and i
    int px, py, pi;
        
    // filter pixel sum
    Vec pixelSum;
    int pixelNum;
        
    // For each pixel in the source
    for (int i = 0; i < pix; ++i) {
        
        // Set Vec to center pixel
        pixelSum = c[i];
        pixelNum = 1;
        
        // Find x and y values
        x = i % w;
        y = ( (h*w) - i - 1) / w;
        
        // Find corner indices
        cx = x - move;
        cy = y - move;
    
        for (int v = 0; v < size; ++v) {
            for (int u = 0; u < size; ++u) {
            
                // Get potential pixel's x and y coords
                px = u + cx;
                py = v + cy;
                
                // Check that index is valid
                if (px >= 0 && px < w && py >= 0 && py < h) {
                    
                    // Compute potential pixel's index
                    pi = (h - py - 1) * w + px;
                    
                    // DEBUG
                    // pixelSum = pixelSum + c[pi] * kernel[v][u];
                    // pixelNum += kernel[v][u];
                    
                    
                    // Check that object id's are the same
                    
                    if (o[pi] == o[i]) {
                        pixelSum = pixelSum + c[pi] * kernel[v][u];
                        pixelNum += kernel[v][u];
                    }
                    
                    
                }           
            }
        } // filterPix has all contributing pixel values
    
        // Compute the median value
        g[i] = pixelSum * (1.0/(double)pixelNum);
    }
    
    // Return the vector of median values
    return g;       

}

// *******************************************************************
// *******************************************************************
// *******************************************************************


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
    int med1 = argc >= 3 ? atoi(argv[2]) : 1;
    int med2 = argc >= 4 ? atoi(argv[3]) : 1;
    int gaus = argc >= 5 ? atoi(argv[4]) : 1;
    // *******************************************************************



    Vec cx = Vec(w*.5135/h), cy = (cx.cross(cam.d)).normalize()*.5135;
    std::vector<Vec> c(w*h);
    
    // *******************************************************************
    std::vector<int> o(w*h);
    // initialize histogram to 0
    int id_count[sizeof(spheres)] = {};
    int id;

#pragma omp parallel for schedule(dynamic, 1)
    for ( int y = 0; y < h; y++ ) {
        for ( int x = 0; x < w; x++ ) {
           
           	// *******************************************************************
           	for (int d = 0; d < sizeof(spheres); ++d)
           		id_count[d] = 0;
           
            const int i = (h - y - 1)*w + x;

            // *******************************************************************
            // COMPUTE OBJECT HERE

            // Shoot ray through the center of pixel i
            double dist;
            id = 0;
            Vec dir = cx*(((x+.5)/w) - .5) + cy*(((y+.5)/h) - .5) + cam.d;
            intersect(Ray(cam.o,dir.normalize()),dist,id);
            o[i] = id;


            for ( int sy = 0; sy < 2; ++sy ) {
                for ( int sx = 0; sx < 2; ++sx ) {
                    Vec r;
                    // *******************************************************************
                    int id_sum = 0;
                    for ( int s = 0; s<samps; s++ ) {
                        double r1 = 2*rng(), dx = r1<1 ? sqrt(r1)-1 : 1-sqrt(2-r1);
                        double r2 = 2*rng(), dy = r2<1 ? sqrt(r2)-1 : 1-sqrt(2-r2);
                        Vec d = cx*(((sx+.5 + dx)/2 + x)/w - .5) + cy*(((sy+.5 + dy)/2 + y)/h - .5) + cam.d;
                        // *******************************************************************
                        r = r + receivedRadiance(Ray(cam.o, d.normalize()), 1, true)*(1./samps);
                        // r = r + receivedRadiance(Ray(cam.o, d.normalize()), 1, true, id)*(1./samps);
                        
                        // double dxy;
                        // intersect(Ray(cam.o, d.normalize()),dxy,id);
                        // id_count[id]++;
                        // std::cout << id << ' ';
                        // o[i] = id;
                        
                    }
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z))*.25;
                }
            }
            
            // *******************************************************************
            
            /*
            o[i] = -1; // std::max_element(id_count,id_count+sizeof(spheres));
            int maxid = 0;
            for (int u = 0; u < sizeof(spheres)-1; ++u) {
                if (id_count[u] > maxid) {
                    o[i] = u;
                    maxid = id_count[u];
                }
            }
            // std::cout << o[i] << ' ';
            */

           // double value = (double)(o[i]+1.0) * (1.0/7.0);
           // Vec newVal(value,value,value);
           // c[i] = newVal;

            
        }
#pragma omp critical
        fprintf(stderr,"\rRendering (%d spp) %6.2f%%",samps*4,100.*y/(h-1));
    }
    fprintf(stderr, "\n");


    // c = boxObjFilter(c, o, h, w, 7);

    // o = medArrFilter(o, h, w, med1);

    c = medObjFilter(c, o, h, w, 7);
    c = gausObjFilter(c, o, h, w, 7);
    c = boxObjFilter(c, o, h, w, 3);

	// Initiate timer
	
	// This will be where mean and median functions will go
	// u[i] = mean(c[i],o[i]), where o[i] is the object of pixel c[i]
	// m[i] = median(c[i],o[i])
	// g[i] = gauss(c[i],o[i])
	
	// After these are implemented test both visual and temporal results
	
	// Next, implement these online, where the algorithm keeps a running track of the current values
	// Test again for time
	
	// receivedRadiance must also return the object number, if no object, then -1
	// this could pass in an array, and index, and receivedRadiance will alter the value inside the
	// method
	
	// Consider making child of Vec which has an extra variable for object number
	
	// For online, need to maintain raw buffer for initial values of c[i]
	
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
