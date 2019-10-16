#pragma once
#include <unordered_set>
#include "automata-base.h"
#include "simplest-cell.h"


namespace cgr {

template <std::size_t Dim, typename Real = default_real>
class simplest_automata 
    : automata_base<Dim, simplest_cell<Dim, Real>, cell_mut_group::mutable_only> {
public:
    using crystallite_type = simplest_cell::crystallite_type;
    using material_type = simplest_crystallite::material_type;
    using orientation_type = crystallite_type::orientation_type;
    using grow_dir = simplest_material::grow_dir;

    bool stop_condition() const override {
        for (auto pcell : *this)
            if (pcell->crystallinity_degre < static_cast<Real>(1.0))
                return false;

        return true;
    }
    bool iterate() override {

    }

    simplest_automata(std::size_t default_range, nbhood_kind default_nbhood_kind)
        : automata_base<Dim, simplest_cell<Dim, Real>, cell_mut_group::mutable_only>(
            default_range, default_nbhood_kind) {}


private:
};

} // namespace cgr
