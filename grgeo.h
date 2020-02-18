// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "geo.h"


namespace cgr::grgeo {

using tag_type = geo::tag_type;
using utag_type = geo::utag_type;
using real_type = geo::real_type;
using vec3r = geo::vec3r;

constexpr std::size_t dim = 3;
using automata_type = cgr::simple_automata<dim, nbhood_kind::euclid, real_type>;
using cell_type = cgr::simple_cell<dim, real_type>;
using cells_container = std::vector<cell_type*>;
using grain_type = cell_type::grain_type;
using grains_container = cell_type::grains_container;

struct gr_point {
    utag_type tag;
    const grains_container* pgrains;

    gr_point(utag_type tag, const grains_container* pgrains)
        : tag(tag), pgrains(pgrains) {}
};

struct gr_line {
    utag_type tag;
    const grains_container* pgrains;

    gr_line(utag_type tag, const grains_container* pgrains)
        : tag(tag), pgrains(pgrains) {}
};

struct gr_surface {
    utag_type tag;
    const grains_container* pgrains;

    gr_surface(utag_type tag, const grains_container* pgrains)
        : tag(tag), pgrains(pgrains) {}
};

struct gr_volume {
    utag_type tag;
    const grain_type* pgrain;

    gr_volume(utag_type tag, const grain_type* pgrain)
        : tag(tag), pgrain(pgrain) {}
};

struct gr_geometry {
    geo::geometry geometry;

    std::vector<gr_volume>  gr_volumes;
    std::vector<gr_surface> gr_surfaces;
    std::vector<gr_line>    gr_lines;
    std::vector<gr_point>   gr_points;

    utag_type add_gr_volume(const grain_type* pgrain, geo::volume::tags_container surface_tags = geo::volume::tags_container()) {
        utag_type tag = geometry.add_volume(std::move(surface_tags));
        gr_volumes.emplace_back(tag, pgrain);
        return tag;
    }
    utag_type add_gr_surface(const grains_container* pgrains, geo::surface::tags_container line_tags = geo::surface::tags_container()) {
        utag_type tag = geometry.add_surface(std::move(line_tags));
        gr_surfaces.emplace_back(tag, pgrains);
        return tag;
    }
    utag_type add_gr_plane_surface(const grains_container* pgrains, geo::surface::tags_container line_tags = geo::surface::tags_container()) {
        utag_type tag = geometry.add_plane_surface(std::move(line_tags));
        gr_surfaces.emplace_back(tag, pgrains);
        return tag;
    }
    utag_type add_gr_line(const grains_container* pgrains, geo::line::tags_container point_tags = geo::line::tags_container()) {
        utag_type tag = geometry.add_line(std::move(point_tags));
        gr_lines.emplace_back(tag, pgrains);
        return tag;
    }
    utag_type add_gr_point(const grains_container* pgrains, vec3r x = vec3r()) {
        utag_type tag = geometry.add_point(x);
        gr_points.emplace_back(tag, pgrains);
        return tag;
    }

    void clear() {
        gr_volumes.clear();
        gr_surfaces.clear();
        gr_lines.clear();
        gr_points.clear();

        gr_volumes.shrink_to_fit();
        gr_surfaces.shrink_to_fit();
        gr_lines.shrink_to_fit();
        gr_points.shrink_to_fit();
    }
};

} // namespace cgr
