#pragma once
#include <vector>
#include <numeric>
#include "simple-grain.h"


namespace cgr {

template <std::size_t Dim, typename Real = double>
struct simple_cell {
    using grain_type = cgr::simple_grain<Dim, Real>;

    Real crystallinity; // 0.0 to 1.0
    std::vector<grain_type*> grains;

    simple_cell(Real crystallinity = 0.0, grain_type* crystallite = nullptr)
        : crystallinity{crystallinity} {
        if (crystallite)
            grains.push_back(crystallite);
    }
};

} // namespace cgr
