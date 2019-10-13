// Copyright © 2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "automata-base.h"


namespace cgr {

template <std::size_t Dim, typename Cell>
class pos_automata_base : automata_base<Dim, Cell> {
public:
    using positons_container = std::unordered_map<Cell*, veci>;

    veci pos(const Cell* cell) const {
        return m_positons[cell];
    }
    void pos(const Cell* cell, const veci& pos) {
        m_positons[cell] = pos;
    }
    virtual void reserve(std::size_t count) {
        m_positons.reserve(count);
    }
    virtual void erase(const Cell* cell) {
        m_positons.erase(cell);
    }

    virtual ~pos_automata_base() {}


private:
    positons_container m_positons;
};

} // namespace cgr
