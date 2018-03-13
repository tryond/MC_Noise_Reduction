#ifndef NORMALFILTER_H
#define NORMALFILTER_H

#include "Vec.h"
#include <vector>
#include <algorithm>
#include <iostream>

std::vector<Vec> normalMean(std::vector<Vec> c, std::vector<int> o, std::vector<Vec> n, int h, int w, int size);

std::vector<std::vector<int>> normalGaussianFilter(int size);

std::vector<Vec> normalGaussian(std::vector<Vec> c, std::vector<int> o, std::vector<Vec> n, int h, int w, int size);

#endif