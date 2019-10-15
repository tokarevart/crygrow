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

template <typename Automata, nbhood_kind NbhoodKind>
class nbh_cont_pos_automata_base : Automata {
public:
    static constexpr auto nbhood_kind = NbhoodKind;
    using nbhood_type = nbhood<NbhoodKind, dim, cell_type>;
    using nbhoods_container = std::unordered_map<
        Cell*, std::unique_ptr<nbhood_type>>;

    std::size_t default_range() const {
        return m_default_range;
    }
    const nbhood* nbhood(const Cell* cell) const {
        auto search = m_nbhoods.find(cell);
        return search != m_nbhoods.end() ? search->second.get() : nullptr;
    }
    void nbhood(const Cell* cell, std::size_t range = default_range()) {
        if (!cell)
            return;

        m_nbhoods[cell] = std::make_unique<nbhood>(
            cell, range,
            [this](const Cell* cell) { return this->pos(cell); },
            [this](const veci& pos) { return this->cell(pos); });
    }
    void erase_nbhood(const Cell* cell) {
        m_nbhoods.erase(cell);
    }
    void reserve_nbhoods(std::size_t count) {
        m_nbhoods.reserve(count);
    }

    nbh_cont_pos_automata_base() = default;
    nbh_cont_pos_automata_base(std::size_t default_range)
        : m_default_range{default_range} {}
    virtual ~nbh_cont_pos_automata_base() {}


private:
    std::size_t m_default_range = 1;
    nbhoods_container m_nbhoods;
};

} // namespace cgr
