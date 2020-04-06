// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "material.h"
#include "cgralgs.h"


namespace cgr {

template <std::size_t Dim, typename Real = double>
class grain {
public:
    using material_type = cgr::material<Dim, Real>;
    using orientation_type = orientation_t<Dim, Real>; // |each vec| == 1

    const material_type* material() const {
        return m_material;
    }
    const orientation_type& orientation() const {
        return m_orientation;
    }

    grain(const material_type* mater, 
          const orientation_type& orien = orientation_type::identity())
        : m_material{ mater }, m_orientation{ orien } {}


private:
    const material_type* m_material;
    orientation_type m_orientation;
};

} // namespace cgr
