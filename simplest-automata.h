#pragma once
#include <unordered_set>
#include "automata-base.h"
#include "simplest-cell.h"


namespace cgr {

template <std::size_t Dim, typename Real = default_real>
class simplest_automata
    : automata_base<Dim, simplest_cell<Dim, Real>, cell_mut_group::mutable_only> {
public:
    using base = automata_base<Dim, simplest_cell<Dim, Real>, cell_mut_group::mutable_only>;
    using cell_type = typename base::cell_type;
    using nbhood_type = typename base::nbhood_type;
    using nbhoods_container = typename base::nbhoods_container;
    using crystallite_type = typename cell_type::crystallite_type;
    using material_type = typename crystallite_type::material_type;
    using orientation_type = typename crystallite_type::orientation_type;
    using grow_dir = typename material_type::grow_dir;

    const nbhood_type* get_direct_nbhood(const cell_type* cell) const {
        auto search = m_direct_nbhoods.find(const_cast<cell_type*>(cell));
        return search != m_direct_nbhoods.end() ? search->second.get() : nullptr;
    }
    void set_direct_nbhood(const cell_type* cell) {
        base::set_nbhood(cell, nbhood_kind::von_neumann, 1);
    }

    bool stop_condition() const override {
        for (auto pcell : *this)
            if (pcell->crystallinity < 1.0)
                return false;

        return true;
    }
    bool iterate() override {
        if (stop_condition())
            return false;

        for (auto pcell : *this) {
            //if (!get_nbhood(pcell))
            //    set_nbhood(pcell);
            if (!get_direct_nbhood(pcell))
                set_direct_nbhood(pcell);
        }

        // test
        for (auto pcell : *this) {
            if (pcell->crystallinity == 1.0)
                continue;

            auto pdir_nbhood = get_direct_nbhood(pcell);
            //auto pcell_pos = pos(pcell);

            for (auto pnb : *pdir_nbhood) {
                if (pnb->crystallinity < 1.0 ||
                    pnb->crystallites.size() != 1)
                    continue;

                if (std::find(pcell->crystallites.begin(),
                              pcell->crystallites.end(),
                              pnb->crystallites.front()) == pcell->crystallites.end())
                    pcell->crystallites.push_back(pnb->crystallites.front());

                pnb->crystallinity += 0.1;
            }
            if (pcell->crystallinity > 1.0)
                pcell->crystallinity = 1.0;
        }

        return true;
    }

    simplest_automata(std::size_t default_range = 1, 
                      nbhood_kind default_nbhood_kind = nbhood_kind::von_neumann)
        : base(default_range, default_nbhood_kind) {}


private:
    nbhoods_container m_direct_nbhoods;
};

} // namespace cgr
