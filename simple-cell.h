#pragma once
#include <vector>
#include <numeric>
#include "simple-grain.h"


namespace cgr {

template <std::size_t Dim, typename Real = double>
struct simple_cell {
    using grain_type = cgr::simple_grain<Dim, Real>;
    using grains_container = std::vector<grain_type*>;

    Real crystallinity; // 0.0 to 1.0
    grains_container grains;

    simple_cell(Real crystallinity = 0.0, grain_type* crystallite = nullptr)
        : crystallinity{crystallinity} {
        if (crystallite)
            grains.push_back(crystallite);
    }
};

} // namespace cgr
