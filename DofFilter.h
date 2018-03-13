#ifndef DOFFILTER_H
#define DOFFILTER_H

#include "Vec.h"
#include <vector>
#include <algorithm>
#include <limits>
#include <stdio.h>
#include <math.h>
#include <iostream>

// levels: amount of levels from which to separate image 
// focLevel: which level should be the center of focus, in [0,levels-1]
// if left undefined by user, function sets focLevel to midpoint of image
std::vector<Vec> dofGray(std::vector<Vec> c, std::vector<int> o, std::vector<int> d, int h, int w, int levels, double focLevel = -1);

std::vector<Vec> dofMean(std::vector<Vec> c, std::vector<int> o, std::vector<int> d, int h, int w, int levels, double focLevel = -1);

#endif