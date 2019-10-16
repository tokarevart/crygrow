#pragma once
#include "automata-base.h"
#include "simplest-cell.h"


namespace cgr {

template <std::size_t Dim, typename Real = default_real>
class simplest_automata 
    : automata_base<Dim, simplest_cell<Dim, Real>, cell_mut_group::mutable_only> {
public:
    using crystallite_type = simplest_cell::crystallite_type;
    using orientation_type = crystallite_type::orientation_type;
    using vec = spt::vec;
    using grow_dir = simplest_material::grow_dir;



private:

};

} // namespace cgr
