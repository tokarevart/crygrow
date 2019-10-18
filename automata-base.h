// Copyright © 2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <memory>
#include <unordered_map>
#include <optional>
#include "neighborhood.h"
#include "cell-mutability.h"
#include "vec.h"


namespace cgr {

template <typename Automata>
class cell_iterator;


template <std::size_t Dim, typename Cell, cell_mut_group CellMutGr>
class automata_base {
public:
    static constexpr std::size_t dim = Dim;
    using cell_type = Cell;
    using veci = spt::veci<Dim>;
    using vecu = spt::vec<Dim, std::uint64_t>;
    using nbhood_type = nbhood<Dim, Cell>;
    // try specify hasher excplicitly if there is an error
    using cells_container = std::unordered_map<veci, std::unique_ptr<Cell>>;
    using positons_container = std::unordered_map<Cell*, veci>;
    using nbhoods_container = std::unordered_map<Cell*, std::unique_ptr<nbhood_type>>;
    using iterator = cell_iterator<automata_base>;
    static constexpr cell_mut_group cell_mut_group = CellMutGr;
    
    iterator begin() const {
        return {m_cells.cbegin()};
    }
    iterator end() const {
        return {m_cells.cend()};
    }

    Cell* cell(const veci& pos) const {
        auto search = m_cells.find(pos);
        return search != m_cells.end() ? search->second.get() : nullptr;
    }
    void  cell(const veci& pos, const Cell* new_cell) {
        if (!new_cell) {
            auto pcell = cell(pos);
            if (pcell)
                erase(pos, pcell);
        } else {
            m_cells[pos].reset(const_cast<Cell*>(new_cell));
            this->pos(new_cell, pos);
        }
    }
    
    veci pos(const Cell* kcell) const {
        return m_positons.at(const_cast<cell_type*>(kcell));
    }
    void pos(const Cell* cell, const veci& pos) {
        m_positons[const_cast<Cell*>(cell)] = pos;
    }

    void reserve(std::size_t count) {
        reserve_cells(count);
        reserve_pos(count);
        reserve_nbhoods(count);
    }
    void reserve_cells(std::size_t count) {
        m_cells.reserve(count);
    }
    void reserve_pos(std::size_t count) {
        m_positons.reserve(count);
    }
    void reserve_nbhoods(std::size_t count) {
        m_nbhoods.reserve(count);
    }

    void erase(const veci& pos, std::optional<const Cell*> corresp_cell) {
        auto pcell = corresp_cell ? corresp_cell.value() : cell(pos);
        erase_cell(pos);
        erase_pos(pcell);
        erase_nbhood(pcell);
    }
    void erase_cell(const veci& pos) {
        m_cells.erase(pos);
    }
    void erase_pos(const Cell* cell) {
        m_positons.erase(const_cast<Cell*>(cell));
    }
    void erase_nbhood(const Cell* cell) {
        m_nbhoods.erase(const_cast<Cell*>(cell));
    }
    
    std::size_t default_range() const {
        return m_default_range;
    }
    std::size_t default_nbhood_kind() const {
        return m_default_nbhood_kind;
    }
    const nbhood_type* get_nbhood(const Cell* cell) const {
        auto search = m_nbhoods.find(cell);
        return search != m_nbhoods.end() ? search->second.get() : nullptr;
    }
    void set_nbhood(const Cell* cell) {
        set_nbhood(cell, default_nbhood_kind(), default_range());
    }
    void set_nbhood(const Cell* cell, std::size_t range) {
        set_nbhood(cell, default_nbhood_kind(), range);
    }
    void set_nbhood(const Cell* cell, nbhood_kind nbhood_kind) {
        set_nbhood(cell, nbhood_kind, default_range());
    }
    void set_nbhood(const Cell* cell, nbhood_kind nbhood_kind, std::size_t range) {
        if (!cell)
            return;

        m_nbhoods[const_cast<cell_type*>(cell)] = std::make_unique<nbhood_type>(
            cell, nbhood_kind, range,
            [this](const cell_type* pcell) { return this->pos(pcell); },
            [this](const veci& pos) { return this->cell(pos); });
    }
    
    virtual bool stop_condition() const = 0;
    virtual bool iterate() = 0;

    automata_base(std::size_t default_range, nbhood_kind default_nbhood_kind)
        : m_default_range{default_range}, m_default_nbhood_kind{default_nbhood_kind} {}

    virtual ~automata_base() {}


private:
    nbhood_kind m_default_nbhood_kind;
    std::size_t m_default_range;

    cells_container m_cells;
    positons_container m_positons;
    nbhoods_container m_nbhoods;
};


template <typename Automata>
struct cell_iterator_base {
    using iterator_category = std::forward_iterator_tag;
    using cell_type = typename Automata::cell_type;
    using from_iterator = typename Automata::cells_container::const_iterator;
};


template <std::size_t Dim, typename Cell, cell_mut_group CellMutGr>
class cell_iterator<automata_base<Dim, Cell, CellMutGr>>
    : cell_iterator_base<automata_base<Dim, Cell, CellMutGr>> {
public:
    using base = cell_iterator_base<automata_base<Dim, Cell, CellMutGr>>;

    typename base::cell_type* to_ptr() const {
        return m_it->second.get();
    }
    cell_iterator next_until(cell_iterator end) {
        return cell_iterator(*this).adv_next_until(end);
    }
    cell_iterator& adv_next_until(cell_iterator end) {
        if (*this == end)
            return end;

        ++m_it;
        if constexpr (std::is_same_v<CellMutGr, cell_mut_group::universal>)
            return to_ptr()->mutability() == cell_mut::constant_cell ? adv_next_until(end) : *this;
        else if constexpr (std::is_same_v<CellMutGr, cell_mut_group::mutable_only>)
            return *this;
        else
            static_assert(false);
    }

    bool operator==(const cell_iterator& other) const {
        return m_it == other.m_it;
    }
    bool operator!=(const cell_iterator& other) const {
        return m_it != other.m_it;
    }
    cell_iterator& operator++() {
        ++m_it;
        if constexpr (CellMutGr == cell_mut_group::universal)
            while(to_ptr()->mutability() == cell_mut::constant_cell)
                ++m_it;
        
        return *this;
    }
    typename base::cell_type* operator*() const {
        return to_ptr();
    }

    cell_iterator(typename base::from_iterator it) : m_it{it} {}


private:
    typename base::from_iterator m_it;
};

} // namespace cgr
