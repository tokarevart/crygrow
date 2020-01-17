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
#include <limits>
#include <algorithm>
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
    const auto& operator[](std::size_t idx) const {
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
    const auto& operator[](std::size_t idx) const {
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
    const auto& operator[](std::size_t idx) const {
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
    template <typename TagType> vec3r make_line_vector(TagType tag) const {
        auto [p0, p1] = get_line_point_tags(tag);
        return get_point(p1).x - get_point(p0).x;
    }
    template <typename TagType> real_type line_length(TagType tag) const {
        return make_line_vector(tag).magnitude();
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

    real_type compute_surface_gmsh_nonplanarity(utag_type tag) const {
        const surface& sur = get_surface(tag);
        std::vector<vec3r> surpoints;
        surpoints.reserve(sur.size());
        for (tag_type ltag : sur)
            surpoints.push_back(get_point(get_line(ltag)[0]).x);

        real_type maxdist2 = 0.0;
        vec3r normal = compute_surface_point_normal(sur[0], sur[1]);
        vec3r mainpoint = surpoints[1];
        vec3r mp_plus_normal = mainpoint + normal;

        for (vec3r point : surpoints) {
            vec3r proj = spt::project_on_normal_line(point, mainpoint, mp_plus_normal);
            real_type dist2 = (proj - mainpoint).magnitude2();
            if (dist2 > maxdist2)
                maxdist2 = dist2;
        }

        return std::sqrt(maxdist2);
    }
    real_type compute_volume_surfaces_gmsh_nonplanarity(utag_type tag) const {
        real_type maxnonpl = 0.0;
        for (tag_type stag : get_volume(tag)) {
            real_type nonpl = compute_surface_gmsh_nonplanarity(std::abs(stag));
            if (nonpl > maxnonpl)
                maxnonpl = nonpl;
        }
        return maxnonpl;
    }
    real_type compute_volume_surfaces_relative_gmsh_nonplanarity(utag_type tag) const {
        real_type minlinelen = std::numeric_limits<real_type>::max();
        for (tag_type stag : get_volume(tag)) {
            for (tag_type ltag : get_surface(stag)) {
                real_type linelen = line_length(ltag);
                if (linelen < minlinelen)
                    minlinelen = linelen;
            }
        }
        return compute_volume_surfaces_gmsh_nonplanarity(tag) / minlinelen;
    }
    real_type compute_gmsh_nonplanarity() const {
        real_type maxnonpl = 0.0;
        for (std::size_t i = 0; i < surfaces.size(); ++i) {
            real_type rnonpl = compute_surface_gmsh_nonplanarity(idx_to_utag(i));
            if (rnonpl > maxnonpl)
                maxnonpl = rnonpl;
        }
        return maxnonpl;
    }
    real_type compute_relative_gmsh_nonplanarity() const {
        real_type maxrnonpl = 0.0;
        for (std::size_t i = 0; i < volumes.size(); ++i) {
            real_type rnonpl = compute_volume_surfaces_relative_gmsh_nonplanarity(idx_to_utag(i));
            if (rnonpl > maxrnonpl)
                maxrnonpl = rnonpl;
        }
        return maxrnonpl;
    }

    real_type compute_surface_worst_nonplanarity(utag_type tag) const {
        const surface& sur = get_surface(tag);
        std::vector<vec3r> surpoints;
        surpoints.reserve(sur.size());
        for (tag_type ltag : sur)
            surpoints.push_back(get_point(get_line(ltag)[0]).x);

        real_type maxdist2 = 0.0;
        vec3r normal = compute_surface_point_normal(sur.back(), sur.front());
        vec3r mainpoint = surpoints.front();
        vec3r mp_plus_normal = mainpoint + normal;
        for (vec3r point : surpoints) {
            vec3r proj = spt::project_on_normal_line(point, mainpoint, mp_plus_normal);
            real_type dist2 = (proj - mainpoint).magnitude2();
            if (dist2 > maxdist2)
                maxdist2 = dist2;
        }
        for (std::size_t i = 0; i < sur.size() - 1; ++i) {
            normal = compute_surface_point_normal(sur[i], sur[i + 1]);
            mainpoint = surpoints[i + 1];
            mp_plus_normal = mainpoint + normal;

            real_type plmaxdist2 = 0.0;
            for (vec3r point : surpoints) {
                vec3r proj = spt::project_on_normal_line(point, mainpoint, mp_plus_normal);
                real_type dist2 = (proj - mainpoint).magnitude2();
                if (dist2 > plmaxdist2)
                    plmaxdist2 = dist2;
            }
            if (plmaxdist2 > maxdist2)
                maxdist2 = plmaxdist2;
        }

        return std::sqrt(maxdist2);
    }
    real_type compute_volume_surfaces_worst_nonplanarity(utag_type tag) const {
        real_type maxnonpl = 0.0;
        for (tag_type stag : get_volume(tag)) {
            real_type nonpl = compute_surface_worst_nonplanarity(std::abs(stag));
            if (nonpl > maxnonpl)
                maxnonpl = nonpl;
        }
        return maxnonpl;
    }
    real_type compute_volume_surfaces_relative_worst_nonplanarity(utag_type tag) const {
        real_type minlinelen = std::numeric_limits<real_type>::max();
        for (tag_type stag : get_volume(tag)) {
            for (tag_type ltag : get_surface(stag)) {
                real_type linelen = line_length(ltag);
                if (linelen < minlinelen)
                    minlinelen = linelen;
            }
        }
        return compute_volume_surfaces_worst_nonplanarity(tag) / minlinelen;
    }
    real_type compute_worst_nonplanarity() const {
        real_type maxnonpl = 0.0;
        for (std::size_t i = 0; i < surfaces.size(); ++i) {
            real_type rnonpl = compute_surface_worst_nonplanarity(idx_to_utag(i));
            if (rnonpl > maxnonpl)
                maxnonpl = rnonpl;
        }
        return maxnonpl;
    }
    real_type compute_relative_worst_nonplanarity() const {
        real_type maxrnonpl = 0.0;
        for (std::size_t i = 0; i < volumes.size(); ++i) {
            real_type rnonpl = compute_volume_surfaces_relative_worst_nonplanarity(idx_to_utag(i));
            if (rnonpl > maxrnonpl)
                maxrnonpl = rnonpl;
        }
        return maxrnonpl;
    }
    
    real_type optimize_surface_nonplanarity(utag_type tag) {
        surface& sur = get_surface(tag);
        std::vector<vec3r> surpoints;
        surpoints.reserve(sur.size());
        for (tag_type ltag : sur)
            surpoints.push_back(get_point(get_line(ltag)[0]).x);

        real_type maxdist2 = 0.0;
        real_type minofmaxdist2 = std::numeric_limits<real_type>::max();
        std::size_t minofmaxdist2idx = 0;
        vec3r normal = compute_surface_point_normal(sur.back(), sur.front());
        vec3r mainpoint = surpoints.front();
        vec3r mp_plus_normal = mainpoint + normal;
        for (vec3r point : surpoints) {
            vec3r proj = spt::project_on_normal_line(point, mainpoint, mp_plus_normal);
            real_type dist2 = (proj - mainpoint).magnitude2();
            if (dist2 > maxdist2)
                maxdist2 = dist2;
            if (dist2 > minofmaxdist2) {
                minofmaxdist2 = dist2;
                minofmaxdist2idx = sur.size() - 1;
            }
        }
        for (std::size_t i = 0; i < sur.size() - 1; ++i) {
            normal = compute_surface_point_normal(sur[i], sur[i + 1]);
            mainpoint = surpoints[i + 1];
            mp_plus_normal = mainpoint + normal;

            real_type plmaxdist2 = 0.0;
            for (vec3r point : surpoints) {
                vec3r proj = spt::project_on_normal_line(point, mainpoint, mp_plus_normal);
                real_type dist2 = (proj - mainpoint).magnitude2();
                if (dist2 > plmaxdist2)
                    plmaxdist2 = dist2;
            }
            if (plmaxdist2 > maxdist2)
                maxdist2 = plmaxdist2;
            if (plmaxdist2 < minofmaxdist2) {
                minofmaxdist2 = plmaxdist2;
                minofmaxdist2idx = i;
            }
        }

        std::rotate(sur.begin(), sur.begin() + minofmaxdist2idx, sur.end());
        return std::sqrt(minofmaxdist2);
    }
    real_type optimize_nonplanarity() {
        real_type maxnonpl = 0.0;
        for (std::size_t i = 0; i < surfaces.size(); ++i) {
            real_type nonpl = optimize_surface_nonplanarity(idx_to_utag(i));
            if (nonpl > maxnonpl)
                maxnonpl = nonpl;
        }
        return maxnonpl;
    }

    real_type compute_surface_best_nonplanarity(utag_type tag) const {
        const surface& sur = get_surface(tag);
        std::vector<vec3r> surpoints;
        surpoints.reserve(sur.size());
        for (tag_type ltag : sur)
            surpoints.push_back(get_point(get_line(ltag)[0]).x);

        real_type minofmaxdist2 = std::numeric_limits<real_type>::max();
        vec3r normal = compute_surface_point_normal(sur.back(), sur.front());
        vec3r mainpoint = surpoints.front();
        vec3r mp_plus_normal = mainpoint + normal;
        for (vec3r point : surpoints) {
            vec3r proj = spt::project_on_normal_line(point, mainpoint, mp_plus_normal);
            real_type dist2 = (proj - mainpoint).magnitude2();
            if (dist2 > minofmaxdist2)
                minofmaxdist2 = dist2;
        }
        for (std::size_t i = 0; i < sur.size() - 1; ++i) {
            normal = compute_surface_point_normal(sur[i], sur[i + 1]);
            mainpoint = surpoints[i + 1];
            mp_plus_normal = mainpoint + normal;

            real_type plmaxdist2 = 0.0;
            for (vec3r point : surpoints) {
                vec3r proj = spt::project_on_normal_line(point, mainpoint, mp_plus_normal);
                real_type dist2 = (proj - mainpoint).magnitude2();
                if (dist2 > plmaxdist2)
                    plmaxdist2 = dist2;
            }
            if (plmaxdist2 < minofmaxdist2)
                minofmaxdist2 = plmaxdist2;
        }

        return std::sqrt(minofmaxdist2);
    }
    real_type compute_volume_surfaces_best_nonplanarity(utag_type tag) const {
        real_type maxnonpl = 0.0;
        for (tag_type stag : get_volume(tag)) {
            real_type nonpl = compute_surface_best_nonplanarity(std::abs(stag));
            if (nonpl > maxnonpl)
                maxnonpl = nonpl;
        }
        return maxnonpl;
    }
    real_type compute_volume_surfaces_relative_best_nonplanarity(utag_type tag) const {
        real_type minlinelen = std::numeric_limits<real_type>::max();
        for (tag_type stag : get_volume(tag)) {
            for (tag_type ltag : get_surface(stag)) {
                real_type linelen = line_length(ltag);
                if (linelen < minlinelen)
                    minlinelen = linelen;
            }
        }
        return compute_volume_surfaces_best_nonplanarity(tag) / minlinelen;
    }
    real_type compute_best_nonplanarity() const {
        real_type maxnonpl = 0.0;
        for (std::size_t i = 0; i < surfaces.size(); ++i) {
            real_type nonpl = compute_surface_best_nonplanarity(idx_to_utag(i));
            if (nonpl > maxnonpl)
                maxnonpl = nonpl;
        }
        return maxnonpl;
    }
    real_type compute_relative_best_nonplanarity() const {
        real_type maxrnonpl = 0.0;
        for (std::size_t i = 0; i < volumes.size(); ++i) {
            real_type rnonpl = compute_volume_surfaces_relative_best_nonplanarity(idx_to_utag(i));
            if (rnonpl > maxrnonpl)
                maxrnonpl = rnonpl;
        }
        return maxrnonpl;
    }

    real_type compute_surface_planarity(utag_type tag) const {
        const surface& sur = get_surface(tag);
        real_type minpl = std::numeric_limits<real_type>::max();

        vec3r rn0 = compute_surface_point_raw_normal(sur.back(), sur.front());
        for (std::size_t i = 0; i < sur.size() - 1; ++i) {
            vec3r rn1 = compute_surface_point_raw_normal(sur[i], sur[i + 1]);
            real_type pl = spt::dot(rn0, rn1) / std::sqrt(rn0.magnitude2() * rn1.magnitude2());
            if (pl < minpl)
                minpl = pl;
        }
        for (std::size_t i = 0; i < sur.size() - 1; ++i) {
            rn0 = compute_surface_point_raw_normal(sur[i], sur[i + 1]);
            for (std::size_t j = 0; j < sur.size() - 1; ++j) {
                vec3r rn1 = compute_surface_point_raw_normal(sur[j], sur[j + 1]);
                real_type pl = spt::dot(rn0, rn1) / std::sqrt(rn0.magnitude2() * rn1.magnitude2());
                if (pl < minpl)
                    minpl = pl;
            }
        }
        return minpl;
    }
    real_type compute_planarity() const {
        real_type minpl = std::numeric_limits<real_type>::max();
        for (std::size_t i = 0; i < surfaces.size(); ++i) {
            real_type pl = compute_surface_planarity(idx_to_utag(i));
            if (pl < minpl)
                minpl = pl;
        }
        return minpl;
    }

    vec3r compute_surface_point_raw_normal(tag_type line0_tag, tag_type line1_tag) const {
        vec3r p1top0 = -make_line_vector(line0_tag);
        vec3r p1top2 = make_line_vector(line1_tag);
        return spt::cross(p1top2, p1top0);
    }
    vec3r compute_surface_point_normal(tag_type line0_tag, tag_type line1_tag) const {
        return compute_surface_point_raw_normal(line0_tag, line1_tag).normalize();
    }
    vec3r compute_plane_surface_raw_normal(tag_type tag) const {
        auto& surface_lines = get_surface(tag).line_tags;
        vec3r rawn = compute_surface_point_raw_normal(surface_lines[0], surface_lines[1]);
        return tag > 0 ? rawn : -rawn;
    }
    vec3r compute_surface_raw_normal(tag_type tag) const {
        const surface& sur = get_surface(tag);
        vec3r rawn = compute_surface_point_raw_normal(sur.back(), sur.front());
        for (std::size_t i = 0; i < sur.size() - 1; ++i)
            rawn += compute_surface_point_raw_normal(sur[i], sur[i + 1]);
        return tag > 0 ? rawn : -rawn;
    }
    vec3r compute_plane_surface_normal(tag_type tag) const {
        return compute_plane_surface_raw_normal(tag).normalize();
    }
    vec3r compute_surface_normal(tag_type tag) const {
        return compute_surface_raw_normal(tag).normalize();
    }
    vec3r compute_auto_surface_normal(tag_type tag) const {
        if (get_surface(tag).is_plane)
            return compute_plane_surface_normal(tag);
        else
            return compute_surface_normal(tag);
    }
    void init_plane_surface_normals() {
        surface_normals.clear();
        surface_normals.reserve(surfaces.size());
        for (std::size_t i = 0; i < surfaces.size(); ++i)
            surface_normals.push_back(compute_plane_surface_normal(idx_to_utag(i)));
    }
    void init_surface_normals() {
        surface_normals.clear();
        surface_normals.reserve(surfaces.size());
        for (std::size_t i = 0; i < surfaces.size(); ++i)
            surface_normals.push_back(compute_surface_normal(idx_to_utag(i)));
    }
    void init_auto_surface_normals() {
        surface_normals.clear();
        surface_normals.reserve(surfaces.size());
        for (std::size_t i = 0; i < surfaces.size(); ++i)
            surface_normals.push_back(compute_auto_surface_normal(idx_to_utag(i)));
    }

    bool surface_contains_line(utag_type sur_tag, utag_type line_tag) const {
        for (auto ltag : get_surface(sur_tag))
            if (std::abs(ltag) == line_tag)
                return true;
        return false;
    }
    utag_type find_2surfaces_common_line(utag_type sur0_tag, utag_type sur1_tag) const {
        const surface& sur0 = get_surface(sur0_tag);
        for (tag_type ltag : sur0)
            if (surface_contains_line(sur1_tag, std::abs(ltag)))
                return std::abs(ltag);
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

    void orient_2surfaces_consistently(tag_type mainsur_tag, tag_type& sur_tag) {
        const surface& mainsur = get_surface(mainsur_tag);
        const surface& sur = get_surface(sur_tag);

        vec3r n0 = get_surface_normal(mainsur_tag);
        vec3r n1 = get_surface_normal(sur_tag);
        vec3r linevec = make_line_vector(find_2surfaces_common_line(std::abs(mainsur_tag), std::abs(sur_tag)));

        vec3r a0 = spt::cross(linevec, n0);
        vec3r a1 = spt::cross(linevec, n1);
        real_type cosa0a1 = spt::cos(a0, a1);

        real_type one_sqrt2 = 1.0 / std::sqrt(2.0);
        if ((cosa0a1 >= -1.0 && cosa0a1 < -one_sqrt2) ||
            (cosa0a1 >= one_sqrt2 && cosa0a1 <= 1)) {
            real_type n0n1 = spt::dot(n0, n1);
            real_type a0a1 = spt::dot(a0, a1);
            if ((n0n1 > 0.0 && a0a1 > 0.0) ||
                (n0n1 < 0.0 && a0a1 < 0.0))
                sur_tag = -sur_tag;
        } else {
            vec3r n0n1 = spt::cross(n0, n1);
            vec3r a0a1 = spt::cross(a0, a1);
            if (spt::dot(n0n1, a0a1) > 0.0)
                sur_tag = -sur_tag;
        }
    }
    void orient_volume_surfaces(utag_type vol_tag) {
        volume& vol = get_volume(vol_tag);
        std::vector<std::size_t> suridxs{ 0 };
        suridxs.reserve(vol.size());
        for (std::size_t i = 0; i < vol.size(); ++i) {
            auto adjsurs = find_adj_surfaces_to_surface_in_volume(vol_tag, std::abs(vol[suridxs[i]]));
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
        return define_entity_str("Surface Loop", vol.tag, vol.surface_tags) + "\n"
            + define_entity_str("Volume", vol.tag, std::vector{ vol.tag });
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
