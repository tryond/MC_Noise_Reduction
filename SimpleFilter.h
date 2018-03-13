#ifndef SIMPLEFILTER_H
#define SIMPLEFILTER_H

#include "Vec.h"
#include <vector>
#include <algorithm>
#include <iostream>

std::vector<Vec> simpleMean(std::vector<Vec> c, int h, int w, int size);

Vec simpleMedianVec(std::vector<Vec> v);

std::vector<Vec> simpleMedian(std::vector<Vec> c, int h, int w, int size);

std::vector<std::vector<int>> simpleGaussianFilter(int size);

std::vector<Vec> simpleGaussian(std::vector<Vec> c, int h, int w, int size);

#endif
