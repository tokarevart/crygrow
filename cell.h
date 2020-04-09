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
    using grains_container = std::vector<const grain_type*>;

    bool crysted;
    grains_container grains;

    cell(const grain_type* grain = nullptr, bool crysted = false)
        : crysted{ crysted } {
        if (grain)
            grains.push_back(grain);
    }
};

} // namespace cgr
