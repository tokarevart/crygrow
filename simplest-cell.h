#pragma once
#include <vector>
#include <numeric>
#include "cell-mutability.h"
#include "simplest-crystallite.h"


namespace cgr {

template <std::size_t Dim, typename Real = default_real>
struct simplest_cell {
    using crystallite_type = simplest_crystallite<Dim, Real>;

    Real crystallinity = 0.0; // 0.0 to 1.0
    std::vector<crystallite_type*> crystallites;
};

} // namespace cgr
