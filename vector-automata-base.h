// Copyright © 2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "nbh-cont-pos-automata-base.h"
#include "vector-no-nbh-automata-base.h"


namespace cgr {

template <std::size_t Dim, typename Cell, cell_mut_group CellMutGr, neighborhood_kind NbhoodKind>
class vector_automata_base : nbh_cont_pos_automata_base<
    vector_no_nbh_automata_base<Dim, Cell, CellMutGr>, NbhoodKind> {
public:

    virtual ~vector_automata_base() {}


private:
};

} // namespace cgr
