#ifndef VEC_H
#define VEC_H

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

struct Vec {
    
	Vec(double x_ = 0, double y_ = 0, double z_ = 0);

	// New Vec Functions ********************************
	
	bool operator< (const Vec &b) const;
	bool operator== (const Vec &b) const;
	Vec& operator= (const Vec &b);
	Vec operator/ (double b) const;
	friend std::ostream& operator<<(std::ostream& out, const Vec &b);
	
	Vec average(const Vec &b) const;
	
	// End New Vec Functions ********************************	 
    
    Vec operator+ (const Vec &b) const;
    Vec operator- (const Vec &b) const;
    Vec operator* (double b) const;
    Vec mult(const Vec &b) const;
    Vec& normalize();
    double dot(const Vec &b) const;
    Vec cross(const Vec &b) const;
  
    double x, y, z;		// x, y, and z components
    double epsilon;		// threshold for equality 
	double mag; 		// magnitude of vector

};

struct Ray {
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}

    // New Ray Functions ********************************

    Ray& operator= (const Ray &b);
    Ray operator+ (const Ray &b) const;
    Ray operator- (const Ray &b) const;

    // End New Vec Functions ********************************	 
};


inline double clamp(double x)   {
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

inline int toInt(double x) {
    return static_cast<int>(std::pow(clamp(x), 1.0/2.2)*255+.5);
}

#endif