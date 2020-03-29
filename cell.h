// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <vector>
#include <numeric>
#include "grain.h"


namespace cgr {

template <std::size_t Dim, typename Real = double>
struct cell {
    using grain_type = cgr::grain<Dim, Real>;
    using grains_container = std::vector<grain_type*>;

    Real crystallinity; // 0.0 to 1.0
    grains_container grains;

    cell(Real crystallinity = 0.0, grain_type* grain = nullptr)
        : crystallinity{crystallinity} {
        if (grain)
            grains.push_back(grain);
    }
};

} // namespace cgr
