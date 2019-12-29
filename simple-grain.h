#pragma once
#include "mat.h"
#include "simple-material.h"


namespace cgr {

template <std::size_t Dim, typename Real = double>
class simple_grain {
public:
    using material_type = cgr::simple_material<Dim, Real>;
    using orientation_type = spt::mat<Dim, Real>; // |each vec| == 1

    const material_type* material() const {
        return m_material;
    }
    const orientation_type& orientation() const {
        return m_orientation;
    }

    simple_grain(const material_type* mater, 
                 const orientation_type& orien = orientation_type::identity())
        : m_material{const_cast<material_type*>(mater)}, m_orientation{orien} {}


private:
    material_type* m_material;
    orientation_type m_orientation;
};

} // namespace cgr
