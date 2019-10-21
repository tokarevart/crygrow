#pragma once
#include <cstddef>
#include <vector>
#include "vec.h"


namespace cgr {

using default_real = double;

template <std::size_t Dim, typename Real = default_real>
class simplest_material {
public:
    using grow_dir = spt::vec<Dim, Real>;

    const std::vector<grow_dir>& grow_dirs() const {
        return m_grow_dirs;
    }

    simplest_material(std::vector<grow_dir> gdirs = {}) {
        if (!gdirs.empty())
            m_grow_dirs = std::move(gdirs);
    }


private:
    std::vector<grow_dir> m_grow_dirs;
};

} // namespace cgr
