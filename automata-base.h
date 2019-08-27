#pragma once
#include <cstddef>
#include <vector>
#include <functional>
#include "vec.h"


namespace cgr {

template<typename Cell>
class cell_iterator : std::iterator<std::forward_iterator_tag, Cell> {
public:
    bool operator!=(cell_iterator const& other) const;
    bool operator==(cell_iterator const& other) const;
    typename cell_iterator::reference operator*() const;
    cell_iterator& operator++();

    cell_iterator(Cell* p);
    // finish later

private:
    Cell* p;
};

// try change template Cell to cell_base class and make benchmarks
template <std::size_t Dim, typename Cell>
class automata_base {
public:
    virtual Cell* get_cell(const spt::vec<Dim, std::size_t>& pos) const = 0;
    virtual void  reset_cell(const spt::vec<Dim, std::size_t>& pos, Cell* ptr = nullptr) = 0;

    virtual cell_iterator begin_iter() const = 0;
    virtual cell_iterator next_iter() const = 0;
    virtual cell_iterator end_iter() const = 0;

    virtual bool stop_condition() const = 0;
    virtual bool iterate() = 0;

    virtual ~automata_base() {}
};

} // namespace cgr
