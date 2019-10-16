#pragma once
#include <vector>
#include "cell-mutability.h"
#include "simplest-crystallite.h"


namespace cgr {

using real_type = double;

struct simplest_cell {
    real_type crystallinity_degree; // 0.0 to 1.0
    real_type crystallization_rate;
    std::vector<simplest_crystallite> contains_crystallites;
    // ...
};

} // namespace cgr
