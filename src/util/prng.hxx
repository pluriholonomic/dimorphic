#pragma once
#include <random>

namespace Dimers {

struct random_int {
    random_int(int n, std::default_random_engine *generator)
        : _n(n), _gen(generator), _dist(0, n - 1) {}

    // For now, we'll make this noncopyable for now
    random_int(const random_int &) = delete;

    int generate() { return _dist(*_gen); }

    int _n;
    std::default_random_engine *_gen;
    std::uniform_int_distribution<int> _dist;
};

struct random_dbl {
    random_dbl(std::default_random_engine *generator)
        : _gen(generator), _dist(0.0, 1.0) {}

    double generate() { return _dist(*_gen); }

    std::default_random_engine *_gen;
    std::uniform_real_distribution<double> _dist;
};

}  // namespace Dimers
