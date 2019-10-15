// Copyright © 2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <numeric>
#include <memory>
#include <optional>
#include "pos-automata-base.h"
#include "cell-mutability.h"
#include "neighborhood.h"


namespace cgr {

template <typename Automata, neighborhood_kind NbhoodKind>
class nbh_cont_pos_automata_base : Automata {
public:
    static constexpr auto neighborhood_kind = NbhoodKind;
    using neighborhood_type = neighborhood<NbhoodKind, dim, cell_type>;
    using neighborhoods_container = std::unordered_map<
        Cell*, std::unique_ptr<neighborhood_type>>;

    std::size_t default_range() const {
        return m_default_range;
    }
    const neighborhood* neighborhood(const Cell* cell) const {
        auto search = m_neighborhoods.find(cell);
        return search != m_neighborhoods.end() ? search->second.get() : nullptr;
    }
    void neighborhood(const Cell* cell, std::size_t range = default_range()) {
        if (!cell)
            return;

        m_neighborhoods[cell] = std::make_unique<neighborhood>(
            cell, range,
            [this](const Cell* cell) { return this->pos(cell); },
            [this](const veci& pos) { return this->cell(pos); });
    }
    void erase_neighborhood(const Cell* cell) {
        m_neighborhoods.erase(cell);
    }
    void reserve_neighborhoods(std::size_t count) {
        m_neighborhoods.reserve(count);
    }

    nbh_cont_pos_automata_base() = default;
    nbh_cont_pos_automata_base(std::size_t default_range)
        : m_default_range{default_range} {}
    virtual ~nbh_cont_pos_automata_base() {}


private:
    std::size_t m_default_range = 1;
    neighborhoods_container m_neighborhoods;
};

} // namespace cgr
