// Copyright © 2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "nbh-cont-pos-automata-base.h"
#include "umap-no-nbh-automata-base.h"


namespace cgr {

template <std::size_t Dim, typename Cell, cell_mut_group CellMutGr, nbhood_kind NbhoodKind>
class umap_automata_base : nbh_cont_pos_automata_base<
    umap_no_nbh_automata_base<Dim, Cell, CellMutGr>, NbhoodKind> {
public:
    void reserve(std::size_t count) {
        reserve_cells(count);
        reserve_pos(count);
        reserve_nbhoods(count);
    }
    
    virtual ~umap_automata_base() {}


private:
};

} // namespace cgr
