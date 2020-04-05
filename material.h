// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <vector>
#include "vec.h"


namespace cgr {

template <std::size_t Dim, typename Real = double>
class material {
public:
    using grow_dir = spt::vec<Dim, Real>;

    const std::vector<grow_dir>& grow_dirs() const {
        return m_grow_dirs;
    }

    material(std::vector<grow_dir> gdirs = {})
        : m_grow_dirs{ std::move(gdirs) } {}


private:
    std::vector<grow_dir> m_grow_dirs;
};

} // namespace cgr
