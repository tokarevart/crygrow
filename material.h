// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <vector>
#include "vec.h"
#include "cgralgs.h"


namespace cgr {

template <std::size_t Dim, typename Real = double>
class material {
public:
    using grow_dir_type = grow_dir_t<Dim, Real>;

    const std::vector<grow_dir_type>& grow_dirs() const {
        return m_grow_dirs;
    }

    material(std::vector<grow_dir_type> gdirs = {})
        : m_grow_dirs{ std::move(gdirs) } {}


private:
    std::vector<grow_dir_type> m_grow_dirs;
};

} // namespace cgr
