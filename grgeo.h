#pragma once
#include "geo.h"


namespace cgr::grgeo {

using tag_type = geo::tag_type;
using utag_type = geo::utag_type;
using real_type = geo::real_type;
using pos_type = geo::pos_type;

constexpr std::size_t dim = 3;
using automata_type = cgr::simple_automata<dim, nbhood_kind::euclid, real_type>;
using cell_type = cgr::simple_cell<dim, real_type>;
using cells_container = std::vector<cell_type*>;
using grain_type = cell_type::grain_type;
using grains_container = cell_type::grains_container;

struct gr_vert {
    utag_type tag;
    const grains_container* pgrains;
};

struct gr_edge {
    tag_type tag;
    const grains_container* pgrains;
};

struct gr_face {
    tag_type tag;
    const grains_container* pgrains;
};

struct gr_volume {
    tag_type tag;
    const grain_type* pgrain;
};

struct gr_geometry {
    geo::geometry geometry;

    std::vector<gr_volume> volumes;
    std::vector<gr_face>   faces;
    std::vector<gr_edge>   edges;
    std::vector<gr_vert>   verts;
};

} // namespace cgr
