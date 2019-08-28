#pragma once
#include <cstddef>
#include <vector>
#include <functional>
#include "vec.h"


namespace cgr {

template <typename Automata>
struct cell_iterator_base : std::iterator<std::forward_iterator_tag, Automata::cell_type> {
    using cell_type = Automata::cell_type;
    using value_type = cell_type*;
    using from_iterator = Automata::cells_container_type::iterator;
};

template <typename Automata>
class cell_iterator;


// try change template Cell to cell_base class and make benchmarks
template <std::size_t Dim, typename Cell>
class automata_base {
public:
    static constexpr std::size_t dim = Dim;
    using cell_type = Cell;
    using veci = spt::veci<Dim>;
    using vecu = spt::vec<Dim, std::uint64_t>;

    virtual Cell* get(const vecu& pos) const = 0;
    virtual void reset(const vecu& pos, Cell* ptr = nullptr) = 0;
    
    virtual bool stop_condition() const = 0;
    virtual bool iterate() = 0;

    virtual ~automata_base() {}
};

} // namespace cgr
