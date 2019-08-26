#pragma once
#include <cstddef>
#include <vector>
#include <functional>
#include "vec.h"


namespace cgr {

template <std::size_t N, typename Cell>
class automata_base {
public:
    virtual Cell* get_cell(const spt::vec<N, std::size_t>& pos) const = 0;
    virtual void reset_cell(const spt::vec<N, std::size_t>& pos, Cell* ptr = nullptr) = 0;

    virtual bool iterate(std::size_t num_iters = 1) final {
        for (auto&& step : m_iteration_steps) 
            step();
    }
    virtual void add_iteration_step(std::function<void()> step) final {
        m_iteration_steps.push_back(step);
    }

    virtual ~automata_base() {}


private:
    std::vector<std::function<void()>> m_iteration_steps;
};

} // namespace cgr
