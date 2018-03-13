#ifndef OBJECTFILTER_H
#define OBJECTFILTER_H

#include "Vec.h"
#include <vector>
#include <algorithm>
#include <iostream>

std::vector<Vec> objectMean(std::vector<Vec> c, std::vector<int> o, int h, int w, int size);

Vec objectMedianVec(std::vector<Vec> v);

std::vector<Vec> objectMedian(std::vector<Vec> c, std::vector<int> o, int h, int w, int size);

std::vector<std::vector<int>> objectGaussianFilter(int size);

std::vector<Vec> objectGaussian(std::vector<Vec> c, std::vector<int> o, int h, int w, int size);

#endif