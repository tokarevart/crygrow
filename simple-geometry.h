#pragma once
#include "simple-automata.h"


namespace cgr {

template <typename Automata>
class simple_geometry {
public:
    static constexpr std::size_t dim = Automata::Dim;
    using cell_type = typename Automata::cell_type;
    using cells_container = typename Automata::cells_container;

    cells_container grain_boundaries(const cells_container& cells) const {
        cells_container res;
        for (auto& pcell : cells)
            if (pcell->grains.size() > 1)
                res.push_back(pcell);
        return res;
    }
    // res_type group_boundaries(const cells_container& boundarycells) const {
    //    
    //}


private:
};

} // namespace cgr
