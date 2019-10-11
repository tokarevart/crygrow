#pragma once
#include "automata-base.h"
#include <array>
#include <numeric>
#include <memory>
#include <unordered_map>


namespace cgr {

template <std::size_t Dim, typename Cell>
class umap_automata_base : automata_base<Dim, Cell> {
public:
    // try specify hasher excplicitly if there is error
    using cells_container_type = std::unordered_map<veci, std::unique_ptr<Cell>>;
    using iterator = cell_iterator<vector_automata_base<Dim, Cell>>;

    iterator begin() const {
        return { m_cells.begin() };
    }
    iterator end() const {
        return { m_cells.end() };
    }

    Cell* get(const veci& pos) const override final {
        auto search = m_cells.find(pos);
        return search != m_cells.end() ? search->second.get() : nullptr;
    }
    void reset(const veci& pos, Cell* ptr = nullptr) override final {
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
class cell_iterator<umap_automata_base<Dim, Cell>> 
    : cell_iterator_base<umap_automata_base<Dim, Cell>> {
public:
    pointer to_ptr() const override final {
        return m_it->second.get();
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
    value_type& operator*() const {
        return *to_ptr();
    }

    cell_iterator(from_iterator it)
        : m_it{ it } {}


private:
    from_iterator m_it;
};

} // namespace cgr
