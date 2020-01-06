#pragma once
#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <array>
#include <vector>
#include <memory>
#include <cmath>
#include "vec.h"
#include "sptops.h"
#include "sptalgs.h"
#include "iteration.h"

namespace geo {

using real_type = double;
using tag_type = std::int_fast32_t;
using utag_type = std::uint_fast32_t;

struct geometry {
    std::vector<volume>  volumes;
    std::vector<surface> surfaces;
    std::vector<line>    lines;
    std::vector<point>   points;

    std::size_t tag_to_idx(tag_type tag) const {
        return std::abs(tag) - 1;
    }

    itr::iteration<tag_type> volume_iteration(tag_type volume_tag) {
        return { 
            0, volumes[tag_to_idx(volume_tag)].size(),
            volume_tag > 0 ? itr::direction::forward : itr::direction::reverse 
        };
    }
    itr::iteration<tag_type> surface_iteration(tag_type surface_tag) {
        return {
            0, surfaces[tag_to_idx(surface_tag)].size(),
            surface_tag > 0 ? itr::direction::forward : itr::direction::reverse
        };
    }
    itr::iteration<utag_type> line_iteration(tag_type line_tag) {
        return {
            0, 2,
            line_tag > 0 ? itr::direction::forward : itr::direction::reverse
        };
    }
};

struct point {
    utag_type tag;
    spt::vec3<real_type> x;
};

struct line {
    using tags_container = std::array<utag_type, 2>;
    tags_container points_tags;
    utag_type tag;

    tags_container::iterator begin() {
        return points_tags.begin();
    }
    tags_container::iterator end() {
        return points_tags.end();
    }
    tags_container::reverse_iterator rbegin() {
        return points_tags.rbegin();
    }
    tags_container::reverse_iterator rend() {
        return points_tags.rend();
    }
    std::size_t size() const {
        return points_tags.size();
    }

    tag_type front() {
        return points_tags.front();
    }
    tag_type back() {
        return points_tags.back();
    }
};

struct surface {
    using tags_container = std::vector<tag_type>;
    tags_container lines_tags;
    utag_type tag;
    bool is_plane = false;

    tags_container::iterator begin() {
        return lines_tags.begin();
    }
    tags_container::iterator end() {
        return lines_tags.end();
    }
    tags_container::reverse_iterator rbegin() {
        return lines_tags.rbegin();
    }
    tags_container::reverse_iterator rend() {
        return lines_tags.rend();
    }
    std::size_t size() const {
        return lines_tags.size();
    }

    tag_type front() {
        return lines_tags.front();
    }
    tag_type back() {
        return lines_tags.back();
    }
};

struct volume {
    using tags_container = std::vector<tag_type>;
    tags_container surfaces_tags;
    utag_type tag;

    tags_container::iterator begin() {
        return surfaces_tags.begin();
    }
    tags_container::iterator end() {
        return surfaces_tags.end();
    }
    tags_container::reverse_iterator rbegin() {
        return surfaces_tags.rbegin();
    }
    tags_container::reverse_iterator rend() {
        return surfaces_tags.rend();
    }
    std::size_t size() const {
        return surfaces_tags.size();
    }

    tag_type front() {
        return surfaces_tags.front();
    }
    tag_type back() {
        return surfaces_tags.back();
    }
};

} // namespace geo
