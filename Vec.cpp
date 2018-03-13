#include "Vec.h"

Vec::Vec(double x_, double y_, double z_) { 
	x = x_; 
	y = y_; 
	z = z_; 
	mag = x_*x_ + y_*y_ + z_*z_; 
	epsilon = 0.00001; 
}


// New Vec Functions ********************************

bool Vec::operator< (const Vec &b) const { 
	return mag < b.mag; 
} 

bool Vec::operator== (const Vec &b) const { 
	return abs(x-b.x) > epsilon ? 0 : 
		   abs(y-b.y) > epsilon ? 0 : 
		   abs(z-b.z) > epsilon ? 0 : 1; 
}

Vec& Vec::operator= (const Vec &b) {
	x = b.x;
	y = b.y;
	z = b.z;
	mag = b.mag;
	return *this;
}

Vec Vec::operator/ (double b) const { 
	return Vec(x/b, y/b, z/b); 
} 

std::ostream& operator<<(std::ostream& out, const Vec &b)
{
    out << "(" << b.x << ", " << b.y << ", " << b.z << ")" << " = " << b.mag;
    return out;
}

Vec Vec::average(const Vec &b) const {
	return Vec((x+b.x)/2.0, (y+b.y)/2.0, (z+b.z)/2.0); 
}	

// End New Vec Functions ********************************


Vec Vec::operator+ (const Vec &b) const { 
	return Vec(x+b.x, y+b.y, z+b.z); 
}

Vec Vec::operator- (const Vec &b) const { 
	return Vec(x-b.x, y-b.y, z-b.z); 
}

Vec Vec::operator* (double b) const { 
	return Vec(x*b, y*b, z*b); 
}

Vec Vec::mult(const Vec &b) const { 
	return Vec(x*b.x, y*b.y, z*b.z); 
}

Vec& Vec::normalize() { 
	return *this = *this * (1.0/std::sqrt(x*x+y*y+z*z)); 
}

double Vec::dot(const Vec &b) const { 
	return x*b.x+y*b.y+z*b.z; 
}

Vec Vec::cross(const Vec &b) const { 
	return Vec(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x); 
}
  

// New Ray Functions ********************************

Ray& Ray::operator= (const Ray &b) {
	o = b.o;
	d = b.d;
	return *this;
}

Ray Ray::operator+ (const Ray &b) const { 
	return Ray(o+b.o, (d+b.d).normalize());
}

Ray Ray::operator- (const Ray &b) const { 
	return Ray(o-b.o, (d-b.d).normalize());
}

// End New Ray Functions ********************************
