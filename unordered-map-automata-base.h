#pragma once
#include "automata-base.h"
#include <array>
#include <numeric>
#include <memory>
#include <unordered_map>


namespace cgr {

template <std::size_t Dim, typename Cell>
class unordered_map_automata_base : virtual automata_base<Dim, Cell> {
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
        return m_cells[pos].get();
    }
    void reset(const veci& pos, Cell* ptr = nullptr) override final {
        m_cells[pos].reset(ptr);
    }

    virtual void reserve(std::size_t count) final {
        m_cells.reserve(count);
    }
    virtual void shrink_to_fit() final {
        for (auto it = m_cells.begin(); it != m_cells.end();) {
            if (!it->second)
                it = c.erase(it);
            else
                ++it;
        }
    }

    virtual ~unordered_map_automata_base() {}


private:
    cells_container_type m_cells;
};


template <std::size_t Dim, typename Cell>
class cell_iterator<unordered_map_automata_base<Dim, Cell>> 
    : cell_iterator_base<unordered_map_automata_base<Dim, Cell>> {
public:
    // implementation

    cell_iterator(from_iterator it)
        : m_it{ it } {}


private:
    from_iterator m_it;
};

} // namespace cgr
