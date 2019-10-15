// Copyright © 2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "umap-nonbh-automata-base.h"
#include "neighborhood.h"


namespace cgr {

template <std::size_t Dim, typename Cell, cell_mut_group CellMutGr, neighborhood_type NbhoodType>
class umap_automata_base : umap_nonbh_automata_base<Dim, Cell, CellMutGr> {
public:
    using neighborhoods_container = std::unordered_map<Cell*, std::unique_ptr<neighborhood<NbhoodType, Dim, Cell>>>;

    const neighborhood* neighborhood(const Cell* cell) const {
        auto search = m_neighborhoods.find(cell);
        return search != m_neighborhoods.end() ? search->second.get() : nullptr;
    }
    void neighborhood(const Cell* cell) {
        neighborhood(cell, m_default_range);
    }
    void neighborhood(const Cell* cell, std::size_t range) {
        if (!cell)
            return;

        if (range == 0)
            m_neighborhoods.erase(cell);
        else
            m_neighborhoods[cell] = std::make_unique<neighborhood>(
                cell, range,
                [this](const Cell* cell) { return this->pos(cell); },
                [this](const veci& pos) { return this->cell(pos); });
    }

    umap_automata_base() = default;
    umap_automata_base(std::size_t default_range) 
        : m_default_range{default_range} {}
    virtual ~umap_automata_base() {}


private:
    std::size_t m_default_range = 1;
    neighborhoods_container m_neighborhoods;
};

} // namespace cgr
