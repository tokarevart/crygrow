#pragma once
#include <vector>
#include "cell-mutability.h"
#include "simplest-crystallite.h"


namespace cgr {

template <std::size_t Dim, typename Real = default_real>
struct simplest_cell {
    using crystallite_type = simplest_crystallite<Dim, Real>;

    // sum is from 0.0 to 1.0
    std::vector<std::pair<crystallite_type*, Real>> crystallinity_degrees;
};

} // namespace cgr
