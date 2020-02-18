// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <vector>
#include "vec.h"


namespace cgr {

enum class material_property {
    isotropic,
    anisotropic
};

template <std::size_t Dim, typename Real = double>
class simple_material {
public:
    using grow_dir = spt::vec<Dim, Real>;

    material_property matproperty() const {
        return m_matproperty;
    }
    const std::vector<grow_dir>& grow_dirs() const {
        return m_grow_dirs;
    }

    simple_material(material_property prop, std::vector<grow_dir> gdirs = {})
        : m_matproperty{prop} {
        if (!gdirs.empty())
            m_grow_dirs = std::move(gdirs);
    }


private:
    material_property m_matproperty;
    std::vector<grow_dir> m_grow_dirs;
};

} // namespace cgr
