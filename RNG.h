#ifndef RNG_H
#define RNG_H

#include <random>
#include <vector>
#include <omp.h>

struct RNG {
    RNG();
    void init(int nworkers); 
    double operator()();
    std::uniform_real_distribution<double> distrb;
    std::vector<std::mt19937> engines;
};

#endif