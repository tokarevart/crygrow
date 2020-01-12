// Copyright © 2020 Artem Tokarev. All rights reserved.
// Licensed under the MIT License.

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
using vec3r = spt::vec3<real_type>;

struct point {
    utag_type tag;
    vec3r x;

    auto begin() {
        return x.x.begin();
    }
    auto begin() const {
        return x.x.begin();
    }
    auto end() {
        return x.x.end();
    }
    auto end() const {
        return x.x.end();
    }
    auto rbegin() {
        return x.x.rbegin();
    }
    auto rbegin() const {
        return x.x.rbegin();
    }
    auto rend() {
        return x.x.rend();
    }
    auto rend() const {
        return x.x.rend();
    }
    std::size_t size() const {
        return x.x.size();
    }

    auto& front() {
        return x.x.front();
    }
    auto& front() const {
        return x.x.front();
    }
    auto& back() {
        return x.x.back();
    }
    auto& back() const {
        return x.x.back();
    }

    auto& operator[](std::size_t idx) {
        return x[idx];
    }
    auto& operator[](std::size_t idx) const {
        return x[idx];
    }

    point(utag_type tag, vec3r x = vec3r())
        : tag(tag), x(x) {}
};

struct line {
    using tags_container = std::array<utag_type, 2>;
    tags_container point_tags{ 0, 0 };
    utag_type tag;

    auto begin() {
        return point_tags.begin();
    }
    auto begin() const {
        return point_tags.begin();
    }
    auto end() {
        return point_tags.end();
    }
    auto end() const {
        return point_tags.end();
    }
    auto rbegin() {
        return point_tags.rbegin();
    }
    auto rbegin() const {
        return point_tags.rbegin();
    }
    auto rend() {
        return point_tags.rend();
    }
    auto rend() const {
        return point_tags.rend();
    }
    std::size_t size() const {
        return point_tags.size();
    }

    auto& front() {
        return point_tags.front();
    }
    auto& front() const {
        return point_tags.front();
    }
    auto& back() {
        return point_tags.back();
    }
    auto& back() const {
        return point_tags.back();
    }

    auto& operator[](std::size_t idx) {
        return point_tags[idx];
    }
    auto& operator[](std::size_t idx) const {
        return point_tags[idx];
    }

    line(utag_type tag, tags_container point_tags = tags_container())
        : tag(tag), point_tags(std::move(point_tags)) {}
};

struct surface {
    using tags_container = std::vector<tag_type>;
    tags_container line_tags;
    utag_type tag;
    bool is_plane = false;

    auto begin() {
        return line_tags.begin();
    }
    auto begin() const {
        return line_tags.begin();
    }
    auto end() {
        return line_tags.end();
    }
    auto end() const {
        return line_tags.end();
    }
    auto rbegin() {
        return line_tags.rbegin();
    }
    auto rbegin() const {
        return line_tags.rbegin();
    }
    auto rend() {
        return line_tags.rend();
    }
    auto rend() const {
        return line_tags.rend();
    }
    std::size_t size() const {
        return line_tags.size();
    }

    auto& front() {
        return line_tags.front();
    }
    auto& front() const {
        return line_tags.front();
    }
    auto& back() {
        return line_tags.back();
    }
    auto& back() const {
        return line_tags.back();
    }

    auto& operator[](std::size_t idx) {
        return line_tags[idx];
    }
    auto& operator[](std::size_t idx) const {
        return line_tags[idx];
    }

    surface(utag_type tag, tags_container line_tags = tags_container())
        : tag(tag), line_tags(std::move(line_tags)) {}
};

struct volume {
    using tags_container = std::vector<tag_type>;
    tags_container surface_tags;
    utag_type tag;

    auto begin() {
        return surface_tags.begin();
    }
    auto begin() const {
        return surface_tags.begin();
    }
    auto end() {
        return surface_tags.end();
    }
    auto end() const {
        return surface_tags.end();
    }
    auto rbegin() {
        return surface_tags.rbegin();
    }
    auto rbegin() const {
        return surface_tags.rbegin();
    }
    auto rend() {
        return surface_tags.rend();
    }
    auto rend() const {
        return surface_tags.rend();
    }
    std::size_t size() const {
        return surface_tags.size();
    }

    auto& front() {
        return surface_tags.front();
    }
    auto& front() const {
        return surface_tags.front();
    }
    auto& back() {
        return surface_tags.back();
    }
    auto& back() const {
        return surface_tags.back();
    }

    auto& operator[](std::size_t idx) {
        return surface_tags[idx];
    }
    auto& operator[](std::size_t idx) const {
        return surface_tags[idx];
    }

    volume(utag_type tag, tags_container surface_tags = tags_container())
        : tag(tag), surface_tags(std::move(surface_tags)) {}
};

struct geometry {
    std::vector<volume>  volumes;
    std::vector<surface> surfaces;
    std::vector<line>    lines;
    std::vector<point>   points;

    std::vector<vec3r>   surface_normals;

    template <typename TagType>
    std::size_t tag_to_idx(TagType tag) const {
        if constexpr (std::is_signed_v<TagType>) {
            return static_cast<std::size_t>(std::abs(tag)) - 1;
        } else {
            return tag - 1;
        }
    }
    utag_type idx_to_utag(std::size_t idx) const {
        return idx + 1;
    }

    utag_type add_volume(volume::tags_container surface_tags = volume::tags_container()) {
        utag_type tag = volumes.size() + 1;
        volumes.emplace_back(tag, std::move(surface_tags));
        return tag;
    }
    utag_type add_surface(surface::tags_container line_tags = surface::tags_container()) {
        utag_type tag = surfaces.size() + 1;
        surfaces.emplace_back(tag, std::move(line_tags));
        return tag;
    }
    utag_type add_plane_surface(surface::tags_container line_tags = surface::tags_container()) {
        utag_type tag = add_surface(std::move(line_tags));
        surfaces[tag_to_idx(tag)].is_plane = true;
        return tag;
    }
    utag_type add_line(line::tags_container point_tags = line::tags_container()) {
        utag_type tag = lines.size() + 1;
        lines.emplace_back(tag, std::move(point_tags));
        return tag;
    }
    utag_type add_point(vec3r x = vec3r()) {
        utag_type tag = points.size() + 1;
        points.emplace_back(tag, x);
        return tag;
    }

    template <typename TagType> volume& get_volume(TagType tag) {
        return volumes[tag_to_idx(tag)];
    }
    template <typename TagType> const volume& get_volume(TagType tag) const {
        return volumes[tag_to_idx(tag)];
    }
    template <typename TagType> surface& get_surface(TagType tag) {
        return surfaces[tag_to_idx(tag)];
    }
    template <typename TagType> const surface& get_surface(TagType tag) const {
        return surfaces[tag_to_idx(tag)];
    }
    template <typename TagType> line& get_line(TagType tag) {
        return lines[tag_to_idx(tag)];
    }
    template <typename TagType> const line& get_line(TagType tag) const {
        return lines[tag_to_idx(tag)];
    }
    template <typename TagType> point& get_point(TagType tag) {
        return points[tag_to_idx(tag)];
    }
    template <typename TagType> const point& get_point(TagType tag) const {
        return points[tag_to_idx(tag)];
    }

    template <typename TagType> vec3r get_surface_normal(TagType tag) const {
        return tag > 0 ? surface_normals[tag_to_idx(tag)] : -surface_normals[tag_to_idx(tag)];
    }

    template <typename TagType> auto get_line_point_tags(TagType tag) const {
        const line& line = get_line(tag);
        return tag > 0 ? std::make_pair(line[0], line[1]) : std::make_pair(line[1], line[0]);
    }

    template <typename TagType> itr::range_iter<std::size_t> volume_iter(TagType tag) {
        return {
            0, get_volume(tag).size(),
            tag > 0 ? itr::dir::forward : itr::dir::reverse
        };
    }
    template <typename TagType> itr::range_iter<std::size_t> surface_iter(TagType tag) {
        return {
            0, get_surface(tag).size(),
            tag > 0 ? itr::dir::forward : itr::dir::reverse
        };
    }
    template <typename TagType> itr::range_iter<std::size_t> line_iter(TagType tag) {
        return {
            0, 2,
            tag > 0 ? itr::dir::forward : itr::dir::reverse
        };
    }

    vec3r make_surface_point_raw_normal(tag_type line0_tag, tag_type line1_tag) const {
        auto [p0, p1_0] = get_line_point_tags(line0_tag);
        auto [p1_1, p2] = get_line_point_tags(line1_tag);
        vec3r p1top0 = get_point(p0).x - get_point(p1_0).x;
        vec3r p1top2 = get_point(p2).x - get_point(p1_1).x;
        return spt::cross(p1top2, p1top0);
    }
    vec3r make_plane_surface_raw_normal(tag_type tag) const {
        auto& surface_lines = get_surface(tag).line_tags;
        vec3r rawn = make_surface_point_raw_normal(surface_lines[0], surface_lines[1]);
        return tag > 0 ? rawn : -rawn;
    }
    vec3r make_surface_raw_normal(tag_type tag) const {
        auto& surface_lines = get_surface(tag).line_tags;
        vec3r rawn = make_surface_point_raw_normal(surface_lines.back(), surface_lines.front());
        for (std::size_t i = 0; i < surface_lines.size() - 1; ++i)
            rawn += make_surface_point_raw_normal(surface_lines[i], surface_lines[i + 1]);
        return tag > 0 ? rawn : -rawn;
    }
    vec3r make_plane_surface_normal(tag_type tag) const {
        return make_surface_raw_normal(tag).normalize();
    }
    vec3r make_surface_normal(tag_type tag) const {
        return make_surface_raw_normal(tag).normalize();
    }
    vec3r make_auto_surface_normal(tag_type tag) const {
        if (get_surface(tag).is_plane)
            return make_plane_surface_normal(tag);
        else
            return make_surface_normal(tag);
    }
    void init_plane_surface_normals() {
        surface_normals.clear();
        surface_normals.reserve(surfaces.size());
        for (std::size_t i = 0; i < surfaces.size(); ++i)
            surface_normals.push_back(make_plane_surface_normal(idx_to_utag(i)));
    }
    void init_surface_normals() {
        surface_normals.clear();
        surface_normals.reserve(surfaces.size());
        for (std::size_t i = 0; i < surfaces.size(); ++i)
            surface_normals.push_back(make_surface_normal(idx_to_utag(i)));
    }
    void init_auto_surface_normals() {
        surface_normals.clear();
        surface_normals.reserve(surfaces.size());
        for (std::size_t i = 0; i < surfaces.size(); ++i)
            surface_normals.push_back(make_auto_surface_normal(idx_to_utag(i)));
    }

    bool surface_contains_line(utag_type sur_tag, utag_type line_tag) const {
        for (auto ltag : get_surface(sur_tag))
            if (std::abs(ltag) == line_tag)
                return true;
        return false;
    }
    auto find_adj_surfaces_to_line_in_volume(utag_type vol_tag, utag_type line_tag) const {
        std::pair<std::size_t, std::size_t> res;
        const volume& vol = get_volume(vol_tag);
        for (std::size_t i = 0; i < vol.size(); ++i) {
            if (surface_contains_line(std::abs(vol[i]), line_tag)) {
                res.first = i;
                break;
            }
        }
        for (std::size_t i = res.first + 1; i < vol.size(); ++i) {
            if (surface_contains_line(std::abs(vol[i]), line_tag)) {
                res.second = i;
                break;
            }
        }
        return res;
    }
    std::vector<std::size_t> find_adj_surfaces_to_surface_in_volume(utag_type vol_tag, utag_type sur_tag) const {
        std::vector<std::size_t> res;
        const volume& vol = get_volume(vol_tag);
        const surface& sur = get_surface(sur_tag);
        for (tag_type ltag : sur) {
            auto [s0idx, s1idx] = find_adj_surfaces_to_line_in_volume(vol_tag, std::abs(ltag));
            if (std::abs(vol[s0idx]) == sur_tag)
                res.push_back(s1idx);
            else 
                res.push_back(s0idx);
        }
        return res;
    }

    void orient_2surfaces_consistently(tag_type mainsur, tag_type& sur) {
        // todo
    }
    void orient_volume_surfaces(utag_type vol_tag) {
        volume& vol = get_volume(vol_tag);
        std::vector<std::size_t> suridxs{ 0 };
        suridxs.reserve(vol.size());
        for (std::size_t i = 0; i < vol.size(); ++i) {
            auto adjsurs = find_adj_surfaces_to_surface_in_volume(vol_tag, vol[suridxs[i]]);
            for (std::size_t adjsuridx : adjsurs) {
                if (std::find(suridxs.begin(), suridxs.end(), adjsuridx) != suridxs.end())
                    continue;

                suridxs.push_back(adjsuridx);
                orient_2surfaces_consistently(vol[suridxs[i]], vol[adjsuridx]);
            }
        }
    }
    void orient_surfaces() {
        for (volume& vol : volumes)
            orient_volume_surfaces(vol.tag);
    }

    void orient_surface_lines(utag_type tag) {
        surface& surface = get_surface(tag);
        for (std::size_t mli = 0; mli < surface.size() - 1; ++mli) {
            auto [mp0, mp1] = get_line_point_tags(surface[mli]);
            for (std::size_t li = mli + 1; li < surface.size(); ++li) {
                auto [p0, p1] = get_line_point_tags(surface[li]);
                if (mp1 != p0 && mp1 != p1)
                    continue;
                if (mp1 == p1)
                    surface[li] = -surface[li];
                std::swap(surface[mli + 1], surface[li]);
                break;
            }
        }
    }
    void orient_lines() {
        for (volume& vol : volumes)
            for (tag_type stag : vol)
                orient_surface_lines(std::abs(stag));
    }

    template <typename ExprType>
    std::string expression_str(ExprType expr) const {
        return std::to_string(expr);
    }
    template <template <typename... Args> typename Container, typename ExprType>
    std::string expression_list_str(const Container<ExprType>& exprs) const {
        std::string res = expression_str(exprs[0]);
        for (std::size_t i = 1; i < exprs.size(); ++i)
            res += ", " + expression_str(exprs[i]);
        return res;
    }
    template <template <typename... Args> typename Container, typename ExprType>
    std::string define_entity_str(std::string name, utag_type tag, const Container<ExprType>& exprs) const {
        return name + "(" + expression_str(tag) + ") = {" + expression_list_str(exprs) + "};";
    }

    std::string point_str(const point& pnt) const {
        return define_entity_str("Point", pnt.tag, std::vector(pnt.begin(), pnt.end()));
    }
    std::string line_str(const line& line) const {
        return define_entity_str("Line", line.tag, std::vector(line.begin(), line.end()));
    }
    std::string surface_str(const surface& sur) const {
        return define_entity_str("Line Loop", sur.tag, sur.line_tags) + "\n" 
            + define_entity_str(sur.is_plane ? "Plane Surface" : "Surface", sur.tag, std::vector{ sur.tag });
    }
    std::string volume_str(const volume& vol) const {
        return define_entity_str("Volume", vol.tag, vol.surface_tags);
    }

    void write_points(std::ostream& os) const {
        for (auto& pnt : points)
            os << point_str(pnt) << std::endl;
    }
    void write_lines(std::ostream& os) const {
        for (auto& line : lines)
            os << line_str(line) << std::endl;
    }
    void write_surfaces(std::ostream& os) const {
        for (auto& sur : surfaces)
            os << surface_str(sur) << std::endl;
    }
    void write_volumes(std::ostream& os) const {
        for (auto& vol : volumes)
            os << volume_str(vol) << std::endl;
    }
    void write(std::ostream& os) const {
        write_points(os);
        write_lines(os);
        write_surfaces(os);
        write_volumes(os);
    }

    void clear() {
        volumes.clear();
        surfaces.clear();
        lines.clear();
        points.clear();

        volumes.shrink_to_fit();
        surfaces.shrink_to_fit();
        lines.shrink_to_fit();
        points.shrink_to_fit();
    }
};

} // namespace geo
