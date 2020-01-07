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
#include "range-iter.h"

namespace geo {

using tag_type = std::int_fast32_t;
using utag_type = std::uint_fast32_t;
using real_type = double;
using pos_type = spt::vec3<real_type>;

struct point {
    utag_type tag;
    pos_type x;

    point(utag_type tag, pos_type x = pos_type())
        : tag(tag), x(x) {}
};

struct line {
    using tags_container = std::array<utag_type, 2>;
    tags_container point_tags;
    utag_type tag;

    tags_container::iterator begin() {
        return point_tags.begin();
    }
    tags_container::iterator end() {
        return point_tags.end();
    }
    tags_container::reverse_iterator rbegin() {
        return point_tags.rbegin();
    }
    tags_container::reverse_iterator rend() {
        return point_tags.rend();
    }
    std::size_t size() const {
        return point_tags.size();
    }

    tag_type front() {
        return point_tags.front();
    }
    tag_type back() {
        return point_tags.back();
    }

    utag_type at(std::size_t idx) const {
        return point_tags[idx];
    }
    utag_type operator[](std::size_t idx) const {
        return at(idx);
    }

    line(utag_type tag, tags_container point_tags = tags_container())
        : tag(tag), point_tags(std::move(point_tags)) {}
};

struct surface {
    using tags_container = std::vector<tag_type>;
    tags_container line_tags;
    utag_type tag;
    bool is_plane = false;

    tags_container::iterator begin() {
        return line_tags.begin();
    }
    tags_container::iterator end() {
        return line_tags.end();
    }
    tags_container::reverse_iterator rbegin() {
        return line_tags.rbegin();
    }
    tags_container::reverse_iterator rend() {
        return line_tags.rend();
    }
    std::size_t size() const {
        return line_tags.size();
    }

    tag_type front() {
        return line_tags.front();
    }
    tag_type back() {
        return line_tags.back();
    }

    tag_type at(std::size_t idx) const {
        return line_tags[idx];
    }
    tag_type operator[](std::size_t idx) const {
        return at(idx);
    }

    surface(utag_type tag, tags_container line_tags = tags_container())
        : tag(tag), line_tags(std::move(line_tags)) {}
};

struct volume {
    using tags_container = std::vector<tag_type>;
    tags_container surface_tags;
    utag_type tag;

    tags_container::iterator begin() {
        return surface_tags.begin();
    }
    tags_container::iterator end() {
        return surface_tags.end();
    }
    tags_container::reverse_iterator rbegin() {
        return surface_tags.rbegin();
    }
    tags_container::reverse_iterator rend() {
        return surface_tags.rend();
    }
    std::size_t size() const {
        return surface_tags.size();
    }

    tag_type front() {
        return surface_tags.front();
    }
    tag_type back() {
        return surface_tags.back();
    }

    tag_type at(std::size_t idx) const {
        return surface_tags[idx];
    }
    tag_type operator[](std::size_t idx) const {
        return at(idx);
    }

    volume(utag_type tag, tags_container surface_tags = tags_container())
        : tag(tag), surface_tags(std::move(surface_tags)) {}
};

struct geometry {
    std::vector<volume>  volumes;
    std::vector<surface> surfaces;
    std::vector<line>    lines;
    std::vector<point>   points;

    template <typename TagType>
    std::size_t tag_to_idx(TagType tag) const {
        if constexpr (std::is_signed_v<TagType>) {
            return std::abs(tag) - 1;
        } else {
            return tag - 1;
        }
    }

    utag_type add_volume(volume::tags_container surface_tags = volume::tags_container()) {
        utag_type tag = points.size() + 1;
        volumes.emplace_back(volumes.size() + 1, std::move(surface_tags));
        return tag;
    }
    utag_type add_surface(surface::tags_container line_tags = surface::tags_container()) {
        utag_type tag = points.size() + 1;
        surfaces.emplace_back(surfaces.size() + 1, std::move(line_tags));
        return tag;
    }
    utag_type add_line(line::tags_container point_tags = line::tags_container()) {
        utag_type tag = points.size() + 1;
        lines.emplace_back(lines.size() + 1, std::move(point_tags));
        return tag;
    }
    utag_type add_point(pos_type x = pos_type()) {
        utag_type tag = points.size() + 1;
        points.emplace_back(tag, x);
        return tag;
    }

    template <typename TagType>
    volume& get_volume(TagType tag) {
        return volumes[tag_to_idx(tag)];
    }
    template <typename TagType>
    surface& get_surface(TagType tag) {
        return surfaces[tag_to_idx(tag)];
    }
    template <typename TagType>
    line& get_line(TagType tag) {
        return lines[tag_to_idx(tag)];
    }
    template <typename TagType>
    point& get_point(TagType tag) {
        return points[tag_to_idx(tag)];
    }

    template <typename TagType>
    itr::range_iter<std::size_t> volume_iteration(TagType tag) {
        return {
            0, get_volume(tag).size(),
            tag > 0 ? itr::dir::forward : itr::dir::reverse
        };
    }
    template <typename TagType>
    itr::range_iter<std::size_t> surface_iteration(TagType tag) {
        return {
            0, get_surface(tag).size(),
            tag > 0 ? itr::dir::forward : itr::dir::reverse
        };
    }
    template <typename TagType>
    itr::range_iter<std::size_t> line_iteration(TagType tag) {
        return {
            0, 2,
            tag > 0 ? itr::dir::forward : itr::dir::reverse
        };
    }
};

} // namespace geo
