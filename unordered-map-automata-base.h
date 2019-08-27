#pragma once
#include "automata-base.h"
#include <array>
#include <numeric>
#include <memory>
#include <unordered_map>


namespace cgr {

template <std::size_t Dim, typename Cell>
class unordered_map_automata_base : automata_base<Dim, Cell> {
public:
    Cell* get_cell(const spt::vec<Dim, std::int64_t>& pos) const override final {
        return m_cells[pos].get();
    }
    void reset_cell(const spt::vec<Dim, std::int64_t>& pos, Cell* ptr = nullptr) override final {
        m_cells[pos].reset(ptr);
    }

    virtual void reserve_cells(std::size_t count) final {
        m_cells.reserve(count);
    }
    virtual void shrink_to_fit_cells() final {
        for (auto it = m_cells.begin(); it != m_cells.end();) {
            if (!it->second)
                it = c.erase(it);
            else
                ++it;
        }
    }

    virtual ~unordered_map_automata_base() {}


private:
    std::unordered_map<
        spt::vec<Dim, std::int64_t>,
        std::unique_ptr<Cell>> m_cells; // try specify hasher excplicitly if there is error
};

} // namespace cgr
