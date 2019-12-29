#pragma once
#include <vector>
#include <numeric>
#include "simple-crystallite.h"


namespace cgr {

template <std::size_t Dim, typename Real = double>
struct simple_cell {
    using crystallite_type = cgr::simple_crystallite<Dim, Real>;

    Real crystallinity; // 0.0 to 1.0
    std::vector<crystallite_type*> grains;

    simple_cell(Real crystallinity = 0.0, crystallite_type* crystallite = nullptr)
        : crystallinity{crystallinity} {
        if (crystallite)
            grains.push_back(crystallite);
    }
};

} // namespace cgr
