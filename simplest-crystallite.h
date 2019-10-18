#pragma once
#include "mat.h"
#include "simplest-material.h"


namespace cgr {

template <std::size_t Dim, typename Real = default_real>
class simplest_crystallite {
public:
    using material_type = simplest_material<Dim, Real>;
    using orientation_type = spt::mat<Dim, Real>; // |each vec| == 1

    const material_type* material() const {
        return m_material;
    }
    const orientation_type& orientation() const {
        return m_orientation;
    }

    simplest_crystallite(const material_type* mater, const orientation_type& orien)
        : m_material{mater}, m_orientation{orien} {}


private:
    material_type* m_material;
    orientation_type m_orientation;
};

} // namespace cgr
