#ifndef BRDFFILTER_H
#define BRDFFILTER_H

#include "Vec.h"
#include <vector>
#include <algorithm>
#include <iostream>

std::vector<Vec> brdfVis(std::vector<Vec> c, std::vector<int> o, std::vector<int> r, int h, int w);

std::vector<Vec> brdfMean(std::vector<Vec> c, std::vector<int> o, std::vector<int> r, int h, int w, int size);

std::vector<std::vector<int>> brdfGaussianFilter(int size);

std::vector<Vec> brdfGaussian(std::vector<Vec> c, std::vector<int> o, std::vector<int> r, int h, int w, int size);

#endif