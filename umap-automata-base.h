// Copyright © 2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <array>
#include <numeric>
#include <memory>
#include <optional>
#include "pos-automata-base.h"
#include "cell-mutability.h"


namespace cgr {

template <std::size_t Dim, typename Cell, cell_mut_group CellMutGr>
class umap_automata_base : pos_automata_base<Dim, Cell> {
public:
    // try specify hasher excplicitly if there is error
    using cells_container = std::unordered_map<veci, std::unique_ptr<Cell>>;
    using positons_container = std::unordered_map<Cell*, veci>;
    using iterator = cell_iterator<umap_automata_base<Dim, Cell, CellMutGr>>;
    static constexpr cell_mut_group cell_mut_group = CellMutGr;

    iterator begin() const {
        return { m_cells.cbegin() };
    }
    iterator end() const {
        return { m_cells.cend() };
    }

    Cell* cell(const veci& pos) const override {
        auto search = m_cells.find(pos);
        return search != m_cells.end() ? search->second.get() : nullptr;
    }
    void  cell(const veci& pos, const Cell* new_cell) override {
        if (!new_cell) {
            auto pcell = cell(pos);
            if (pcell)
                erase(pos, pcell);
        } else {
            m_cells[pos].reset(new_cell);
            this->pos(new_cell, pos);
        }
    }
    
    void reserve(std::size_t count) override {
        m_cells.reserve(count);
        pos_automata_base<Dim, Cell>::reserve(count);
    }

    virtual ~umap_automata_base() {}


private:
    cells_container m_cells;

    void erase(const veci& pos, std::optional<const Cell*> corresp_cell) override {
        auto pcell = corresp_cell ? corresp_cell.value() : cell(pos);
        m_cells.erase(pos);
        pos_automata_base<Dim, Cell>::erase(pcell);
    }
};


template <std::size_t Dim, typename Cell>
class cell_iterator<umap_automata_base<Dim, Cell, cell_mut_group::mutable_only>>
    : cell_iterator_base<umap_automata_base<Dim, Cell, cell_mut_group::mutable_only>> {
public:
    cell_type* to_ptr() const override {
        return m_it->second.get();
    }
    cell_iterator next_until(cell_iterator end) {
        return cell_iterator(*this).adv_next_until(end);
    }
    cell_iterator& adv_next_until(cell_iterator end) {
        if (*this == end)
            return end;

        ++m_it;
        return *this;
    }

    bool operator==(const cell_iterator& other) const {
        return m_it == other.m_it;
    }
    bool operator!=(const cell_iterator& other) const {
        return m_it != other.m_it;
    }
    cell_iterator& operator++() const {
        ++m_it;
        return *this;
    }
    cell_type& operator*() const {
        return *to_ptr();
    }

    cell_iterator(from_iterator it) : m_it{it} {}


private:
    from_iterator m_it;
};


template <std::size_t Dim, typename Cell>
class cell_iterator<umap_automata_base<Dim, Cell, cell_mut_group::universal>>
    : cell_iterator_base<umap_automata_base<Dim, Cell, cell_mut_group::universal>> {
public:
    cell_type* to_ptr() const override {
        return m_it->second.get();
    }
    cell_iterator next_until(cell_iterator end) {
        return cell_iterator(*this).adv_next_until(end);
    }
    cell_iterator& adv_next_until(cell_iterator end) {
        if (*this == end)
            return end;

        ++m_it;
        return to_ptr()->mutability() == cell_mut::constant_cell ? adv_next_until(end) : *this;
    }

    bool operator==(const cell_iterator& other) const {
        return m_it == other.m_it;
    }
    bool operator!=(const cell_iterator& other) const {
        return m_it != other.m_it;
    }
    cell_type& operator*() const {
        return *to_ptr();
    }

    cell_iterator(from_iterator it) : m_it{it} {}


private:
    from_iterator m_it;
};

} // namespace cgr
