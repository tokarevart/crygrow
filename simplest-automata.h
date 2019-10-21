#pragma once
#include <unordered_set>
#include <algorithm>
#include <execution>
#include "automata-base.h"
#include "simplest-cell.h"


namespace cgr {

template <std::size_t Dim, typename Real = default_real>
class simplest_automata
    : public automata_base<Dim, simplest_cell<Dim, Real>, cell_mut_group::mutable_only> {
public:
    static constexpr Real epsilon = std::numeric_limits<Real>::epsilon();
    using base = automata_base<Dim, simplest_cell<Dim, Real>, cell_mut_group::mutable_only>;
    using veci = typename base::veci;
    using cell_type = typename base::cell_type;
    using nbhood_type = typename base::nbhood_type;
    using nbhood_pos_type = typename base::nbhood_pos_type;
    using nbhoods_pos_container = typename base::nbhoods_pos_container;
    using cells_delta_container = std::vector<Real>;
    using crystallite_type = typename cell_type::crystallite_type;
    using material_type = typename crystallite_type::material_type;
    using orientation_type = typename crystallite_type::orientation_type;
    using grow_dir = typename material_type::grow_dir;

    bool stop_condition() const override {
        for (auto pcell : *this)
            if (pcell->crystallinity < 1.0)
                return false;

        return true;
    }
    bool iterate() override {
        if (stop_condition())
            return false;

        for (auto& delta : m_cells_delta)
            delta = 0.0;
        
        for (std::size_t i = 0; i < num_cells(); ++i) {
            auto curpos = get_pos(i);
            auto pcell = get_cell(i);
            if (std::abs(pcell->crystallinity - 1.0) <= epsilon * (pcell->crystallinity + 1))
                continue;
            
            auto nbhood = get_nbhood(curpos);
            std::size_t num_acc = 0;
            for (auto pnb : nbhood) {
                if (pnb->crystallinity < 1.0 - epsilon ||
                    pnb->crystallites.size() != 1)
                    continue;

                if (std::find(pcell->crystallites.begin(),
                              pcell->crystallites.end(),
                              pnb->crystallites.front()) == pcell->crystallites.end())
                    pcell->crystallites.push_back(pnb->crystallites.front());

                num_acc++;
            }
            if (num_acc > 0)
                m_cells_delta[i] += static_cast<Real>(num_acc) / 20;
        }
        for (std::size_t i = 0; i < num_cells(); ++i) {
            auto pcell = get_cell(i);
            pcell->crystallinity += m_cells_delta[i];
            m_cells_delta[i] = 0.0;
            if (pcell->crystallinity > 1.0 + epsilon)
                pcell->crystallinity = 1.0;
        }

        return true;
    }

    simplest_automata(const veci& corner0, const veci& corner1, 
                      std::size_t default_range = 1,
                      nbhood_kind default_nbhood_kind = nbhood_kind::euclid)
        : base(corner0, corner1, default_range, default_nbhood_kind) {
        m_cells_delta.assign(m_cells.size(), 0.0);
    }


private:
    cells_delta_container m_cells_delta;

    void initialize_nbhoods_pos() {
        // set_nbhood_pos
    }
};

} // namespace cgr
