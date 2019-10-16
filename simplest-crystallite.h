#pragma once
#include "mat.h"
#include "simplest-material.h"


namespace cgr {

template <std::size_t Dim, typename Real = default_real>
class simplest_crystallite {
public:
    // TODO: use spt::mat instead of std::array
    using orientation_type = spt::mat<Dim, Real>; // |each vec| == 1

    const simplest_material* material() const {
        return m_material;
    }
    const orientation_type& orientation() const {
        return m_orientation;
    }

    simplest_crystallite(const simplest_material* mater, const orientation_type& orien)
        : m_material{mater}, m_orientation{orien} {}


private:
    simplest_material* m_material;
    orientation_type m_orientation;
};

} // namespace cgr
