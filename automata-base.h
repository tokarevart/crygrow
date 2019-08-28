#pragma once
#include <cstddef>
#include <vector>
#include <functional>
#include "vec.h"


namespace cgr {

template<typename Cell>
using cell_iterator_base = std::iterator<std::forward_iterator_tag, Cell>;

template<typename Cell>
class cell_iterator;

// try change template Cell to cell_base class and make benchmarks
template <std::size_t Dim, typename Cell>
class automata_base {
public:
    virtual Cell* get(const spt::vec<Dim, std::size_t>& pos) const = 0;
    virtual void  reset(const spt::vec<Dim, std::size_t>& pos, Cell* ptr = nullptr) = 0;

    virtual cell_iterator<Cell> begin() const = 0;
    virtual cell_iterator<Cell> end() const = 0;

    virtual bool stop_condition() const = 0;
    virtual bool iterate() = 0;

    virtual ~automata_base() {}
};

} // namespace cgr
