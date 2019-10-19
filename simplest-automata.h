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
    using base = automata_base<Dim, simplest_cell<Dim, Real>, cell_mut_group::mutable_only>;
    using veci = typename base::veci;
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
        if (!cell)
            return;

        m_direct_nbhoods[const_cast<cell_type*>(cell)] = std::make_unique<nbhood_type>(
            cell, nbhood_kind::von_neumann, 1,
            [this](const cell_type* pcell) { return this->pos(pcell); },
            [this](const veci& pos) { return this->cell(pos); });
    }

    void reserve(std::size_t count) {
        base::reserve(count);
        reserve_direct_nbhoods(count);
    }
    void reserve_direct_nbhoods(std::size_t count) {
        m_direct_nbhoods.reserve(count);
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
            if (!this->get_nbhood(pcell))
                this->set_nbhood(pcell);
        }
        
        std::for_each(std::execution::par, this->raw_begin(), this->raw_end(), 
        [this](const std::pair<const veci, std::unique_ptr<cell_type>>& l_pair) mutable -> void {
        //for (auto pcell : *this) {
            cell_type* pcell = l_pair.second.get();
            if (std::abs(pcell->crystallinity - 1.0) <= std::numeric_limits<Real>::epsilon() * (pcell->crystallinity + 1))
                //continue;
                return;
            
            auto pnbhood = this->get_nbhood(pcell);
            std::size_t num_acc = 0;
            for (auto pnb : *pnbhood) {
                if (pnb->crystallinity < 1.0 ||
                    pnb->crystallites.size() != 1)
                    continue;

                if (std::find(pcell->crystallites.begin(),
                              pcell->crystallites.end(),
                              pnb->crystallites.front()) == pcell->crystallites.end())
                    pcell->crystallites.push_back(pnb->crystallites.front());

                num_acc++;
            }
            if (num_acc > 0)
                pcell->crystallinity += 0.005 * num_acc;
            if (pcell->crystallinity > 1.0)
                pcell->crystallinity = 1.0;
        });

        return true;
    }

    simplest_automata(std::size_t default_range = 1, 
                      nbhood_kind default_nbhood_kind = nbhood_kind::von_neumann)
        : base(default_range, default_nbhood_kind) {}


private:
    nbhoods_container m_direct_nbhoods;
};

} // namespace cgr
