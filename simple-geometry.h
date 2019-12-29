#pragma once
#include "simple-automata.h"


namespace cgr {

template <nbhood_kind NbhoodKind = nbhood_kind::euclid, typename Real = double>
class simple_geometry {
public:
    static constexpr std::size_t dim = 3;
    using automata_type = cgr::simple_automata<dim, NbhoodKind, Real>;
    using cell_type = typename automata_type::cell_type;
    using cells_container = typename automata_type::cells_container;

    cells_container grain_boundaries(const cells_container& cells) const {
        cells_container res;
        for (auto& pcell : cells)
            if (pcell->grains.size() > 1)
                res.push_back(pcell);
        return res;
    }
    std::array<cells_container, dim> group_boundaries(const cells_container& boundarycells) const {
        
    }


private:
};

} // namespace cgr
