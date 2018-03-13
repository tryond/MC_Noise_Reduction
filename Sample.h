#ifndef SAMPLE_H
#define SAMPLE_H

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "Vec.h"
#include "RNG.h"
#include "Sphere.h"

const double PI = 3.1415926535897932384626433832795;

Vec uniformRandomPSA(const Vec &n);

void luminaireSample(const Sphere &lum, Vec &y, Vec &n, double &pdf);

#endif
