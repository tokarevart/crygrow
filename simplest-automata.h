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
    using nbhoods_container = typename base::nbhoods_container;
    using cells_delta_container = std::unordered_map<cell_type*, Real>;
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

        for (auto pcell : *this) {
            m_cells_delta[pcell] = 0.0;
            if (!this->get_nbhood(pcell))
                this->set_nbhood(pcell);
        }
        
        std::for_each(std::execution::par, this->raw_begin(), this->raw_end(), 
        [this](const std::pair<const veci, std::unique_ptr<cell_type>>& l_pair) mutable -> void {
        //for (auto pcell : *this) {
            cell_type* pcell = l_pair.second.get();
            if (std::abs(pcell->crystallinity - 1.0) <= epsilon * (pcell->crystallinity + 1))
                //continue;
                return;
            
            auto pnbhood = this->get_nbhood(pcell);
            std::size_t num_acc = 0;
            for (auto pnb : *pnbhood) {
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
                m_cells_delta[pcell] += static_cast<Real>(num_acc) / (2 * pnbhood->size());
        });
        for (auto pcell : *this) {
            pcell->crystallinity += m_cells_delta[pcell];
            m_cells_delta[pcell] = 0.0;
            if (pcell->crystallinity > 1.0 + epsilon)
                pcell->crystallinity = 1.0;
        }

        return true;
    }

    simplest_automata(std::size_t default_range = 1, 
                      nbhood_kind default_nbhood_kind = nbhood_kind::von_neumann)
        : base(default_range, default_nbhood_kind) {}


private:
    cells_delta_container m_cells_delta;
};

} // namespace cgr
