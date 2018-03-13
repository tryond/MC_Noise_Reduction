#include "RNG.h"

RNG::RNG() : distrb(0.0, 1.0), engines() {}

void RNG::init(int nworkers) {
    std::random_device rd;
    engines.resize(nworkers);
    for ( int i = 0; i < nworkers; ++i )
        engines[i].seed(rd());
}

double RNG::operator()() {
    int id = omp_get_thread_num();
    return distrb(engines[id]);
}