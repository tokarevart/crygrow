#pragma once
#include <cstddef>
#include <vector>
#include <functional>
#include "vec.h"


namespace cgr {

template <std::size_t Dim, typename Cell>
class automata_base {
public:
    virtual Cell* get_cell(const spt::vec<Dim, std::size_t>& pos) const = 0;
    virtual void reset_cell(const spt::vec<Dim, std::size_t>& pos, Cell* ptr = nullptr) = 0;

    virtual bool stop_condition() const = 0;
    virtual bool iterate() final {
        for (auto& step : m_steps) 
            step();
        return !stop_condition();
    }
    virtual void add_step(std::function<bool()> step) final {
        m_steps.push_back(step);
    }

    virtual ~automata_base() {}


private:
    std::vector<std::function<bool()>> m_steps;
};

} // namespace cgr
