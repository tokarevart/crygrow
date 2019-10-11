#pragma once
#include <array>
#include <unordered_map>
#include <numeric>
#include <memory>
#include "automata-base.h"
#include "cell-mutability.h"


namespace cgr {

template <std::size_t Dim, typename Cell, cell_mut_group CellMutGr>
class umap_automata_base : automata_base<Dim, Cell> {
public:
    // try specify hasher excplicitly if there is error
    using cells_container_type = std::unordered_map<veci, std::unique_ptr<Cell>>;
    using iterator = cell_iterator<umap_automata_base<Dim, Cell, CellMutGr>>;
    static constexpr cell_mut_group cell_mut_group = CellMutGr;

    iterator begin() const {
        return { m_cells.begin() };
    }
    iterator end() const {
        return { m_cells.end() };
    }

    Cell* get(const veci& pos) const override {
        auto search = m_cells.find(pos);
        return search != m_cells.end() ? search->second.get() : nullptr;
    }
    void reset(const veci& pos, Cell* ptr = nullptr) override {
        m_cells[pos].reset(ptr);
    }

    void reserve(std::size_t count) {
        m_cells.reserve(count);
    }
    void shrink_to_fit() {
        for (auto it = m_cells.begin(); it != m_cells.end();) {
            if (!it->second)
                it = c.erase(it);
            else
                ++it;
        }
    }

    virtual ~umap_automata_base() {}


private:
    cells_container_type m_cells;
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

    cell_iterator(from_iterator it) : m_it{ it } {}


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

    cell_iterator(from_iterator it) : m_it{ it } {}


private:
    from_iterator m_it;
};

} // namespace cgr
