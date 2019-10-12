#pragma once
#include <cstddef>
#include <unordered_map>
#include "vec.h"


namespace cgr {

template <typename Automata>
struct cell_iterator_base {
    using iterator_category = std::forward_iterator_tag;
    using cell_type = Automata::cell_type;
    using from_iterator = Automata::cells_container_type::iterator;

    virtual cell_type* to_ptr() const = 0;

    virtual ~cell_iterator_base() {}
};

template <typename Automata>
class cell_iterator;


template <std::size_t Dim, typename Cell>
class automata_base {
public:
    static constexpr std::size_t dim = Dim;
    using cell_type = Cell;
    using veci = spt::veci<Dim>;
    using vecu = spt::vec<Dim, std::uint64_t>;

    virtual Cell* cell(const veci& pos) const = 0;
    virtual void  cell(const veci& pos, const Cell* new_cell) = 0;
    
    virtual bool stop_condition() const = 0;
    virtual bool iterate() = 0;

    virtual ~automata_base() {}
};

} // namespace cgr
