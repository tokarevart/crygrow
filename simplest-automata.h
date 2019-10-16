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

    const nbhood* get_direct_nbhood(const Cell* cell) const {
        auto search = m_direct_nbhoods.find(cell);
        return search != m_direct_nbhoods.end() ? search->second.get() : nullptr;
    }
    void set_direct_nbhood(const Cell* cell) {
        set_nbhood(cell, nbhood_kind::von_neumann, 1);
    }

    bool stop_condition() const override {
        for (auto pcell : *this)
            if (pcell->crystallinity_degre < static_cast<Real>(1.0))
                return false;

        return true;
    }
    bool iterate() override {
        if (stop_condition())
            return false;

        for (auto pcell : *this) {
            set_nbhood(pcell);
            set_direct_nbhood(pcell);
        }

        for (auto pcell : *this) {
            auto pdir_nbhood = get_direct_nbhood(pcell);
            auto pcell_pos = pos(pcell);
            // ...
        }
    }

    simplest_automata(std::size_t default_range, nbhood_kind default_nbhood_kind)
        : automata_base<Dim, simplest_cell<Dim, Real>, cell_mut_group::mutable_only>(
            default_range, default_nbhood_kind) {}


private:
    nbhoods_container m_direct_nbhoods;
};

} // namespace cgr
