// Copyright © 2018-2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

// With help of David Eberly's and Dan Sunday's works.

#pragma once
#include "sptops.h"

// todo: try pass by classes-containers instead of reference and make benchmark
// todo: try vectorization (SIMD) and make benchmark
namespace spt {

template <typename Real>
mat3<Real> rotation(const vec3<Real>& v, Real angle) {
    const Real si = std::sin(angle);
    const Real co = std::cos(angle);
    const Real ic = 1.0 - co;

    return {
        v[0] * v[0] * ic + co, v[1] * v[0] * ic - si * v[2], v[2] * v[0] * ic + si * v[1],
        v[0] * v[1] * ic + si * v[2], v[1] * v[1] * ic + co, v[2] * v[1] * ic - si * v[0],
        v[0] * v[2] * ic - si * v[1], v[1] * v[2] * ic + si * v[0], v[2] * v[2] * ic + co };
}

template <typename ValueType>
bool weak_between(ValueType boundary0, ValueType boundary1, ValueType value) {
    auto l_eps = std::numeric_limits<ValueType>::epsilon() * std::abs(value);
    return
        (value >= std::min(boundary0, boundary1) - l_eps) &&
        (value <= std::max(boundary0, boundary1) + l_eps);
}

template <std::size_t Dim, typename ValueType>
bool weak_in_cuboid(const vec<Dim, ValueType>& corner0,
                    const vec<Dim, ValueType>& corner1,
                    const vec<Dim, ValueType>& point) {
    for (std::size_t i = 0; i < Dim; ++i)
        if (!weak_between(corner0.x[i], corner1.x[i], point.x[i]))
            return false;

    return true;
}

// NOTE: using with vec of integers needs different implementation
template <std::size_t Dim, typename Real>
vec<Dim, Real> project_on_vec(const vec<Dim, Real>& v, const vec<Dim, Real>& on_v) {
    return on_v * (dot(v, on_v) / on_v.magnitude2());
}

template <std::size_t Dim, typename Real>
vec<Dim, Real> project_on_line(
    const vec<Dim, Real>& point,
    const vec<Dim, Real>& line_p0, const vec<Dim, Real>& line_p1) {

    return line_p0 + project_on_vec(point - line_p0, line_p1 - line_p0);
}

// todo: return std::optional<vec<Dim, Real>>
template <std::size_t Dim, typename Real>
bool project_on_segm(
    vec<Dim, Real>& out,
    const vec<Dim, Real>& point,
    const vec<Dim, Real>& segm_p0, const vec<Dim, Real>& segm_p1) {

    auto res = segm_p0 + project_on_vec(point - segm_p0, segm_p1 - segm_p0);

    if (weak_in_cuboid(segm_p0, segm_p1, res)) {
        out = res;
        return true;
    }

    return false;
}

template <typename Real>
vec<3, Real> project_on_plane(
    const vec<3, Real>& v, 
    const vec<3, Real>& plane_v0, const vec<3, Real>& plane_v1) {
    return v - project_on_vec(v, cross(plane_v0, plane_v1));
}

template <typename Real>
vec<3, Real> project_on_plane(
    const vec<3, Real>& point,
    const vec<3, Real>& plane_p0, const vec<3, Real>& plane_p1, const vec<3, Real>& plane_p2) {

    return plane_p0 + project_on_plane(point - plane_p0, plane_p1 - plane_p0, plane_p2 - plane_p0);
}

template <typename Real>
bool does_ray_intersect_plane(
    const vec<3, Real>& dir,
    const vec<3, Real>& pl_p0, const vec<3, Real>& pl_p1, const vec<3, Real>& pl_p2) {

    std::array edges{ pl_p1 - pl_p0, pl_p2 - pl_p0 };

    auto pvec = spt::cross(dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    return std::abs(det) > std::numeric_limits<Real>::epsilon();
}

// todo: return std::optional<vec<3, Real>>
template <typename Real>
bool ray_intersect_plane(
    vec<3, Real>& out_intersect_point,
    const vec<3, Real>& origin, const vec<3, Real>& dir,
    const vec<3, Real>& pl_p0, const vec<3, Real>& pl_p1, const vec<3, Real>& pl_p2) {

    std::array edges{ pl_p1 - pl_p0, pl_p2 - pl_p0 };

    auto pvec = spt::cross(dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    if (std::abs(det) <= std::numeric_limits<Real>::epsilon())
        return false;

    auto tvec = origin - pl_p0;
    auto qvec = spt::cross(tvec, edges[0]);

    auto t = spt::dot(edges[1], qvec) / det;
    out_intersect_point = origin + dir * t;
    return t > static_cast<Real>(0);
}

template <typename Real>
bool does_ray_intersect_triangle(
    const vec<3, Real>& origin, const vec<3, Real>& dir,
    const vec<3, Real>& tr_p0, const vec<3, Real>& tr_p1, const vec<3, Real>& tr_p2) {

    std::array edges{ tr_p1 - tr_p0, tr_p2 - tr_p0 };

    auto pvec = spt::cross(dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    if (std::abs(det) <= std::numeric_limits<Real>::epsilon())
        return false;

    auto inv_det = static_cast<Real>(1) / det;

    auto tvec = origin - tr_p0;
    auto u = spt::dot(tvec, pvec) * inv_det;
    if (u < static_cast<Real>(0) || u > static_cast<Real>(1))
        return false;

    auto qvec = spt::cross(tvec, edges[0]);
    auto v = spt::dot(dir, qvec) * inv_det;
    if (v < static_cast<Real>(0) || u + v > static_cast<Real>(1))
        return false;

    auto t = spt::dot(edges[1], qvec) * inv_det;
    return t >= static_cast<Real>(0);
}

// todo: return std::optional<vec<3>>
template <typename Real>
bool line_intersect_plane(
    vec<3, Real>& out_intersect_point,
    const vec<3, Real>& line_point, const vec<3, Real>& line_dir,
    const vec<3, Real>& plane_p0, const vec<3, Real>& plane_p1, const vec<3, Real>& plane_p2) {

    std::array edges{ plane_p1 - plane_p0, plane_p2 - plane_p0 };

    auto pvec = spt::cross(line_dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    if (std::abs(det) <= std::numeric_limits<Real>::epsilon())
        return false;

    auto tvec = line_point - plane_p0;
    auto qvec = spt::cross(tvec, edges[0]);

    auto t = spt::dot(edges[1], qvec) / det;
    out_intersect_point = line_point + line_dir * t;
    return true;
}

template <typename Real>
vec<3, Real> line_intersect_plane(
    const vec<3, Real>& line_point, const vec<3, Real>& line_dir,
    const vec<3, Real>& plane_p0, const vec<3, Real>& plane_p1, const vec<3, Real>& plane_p2) {

    std::array edges{ plane_p1 - plane_p0, plane_p2 - plane_p0 };

    auto pvec = spt::cross(line_dir, edges[1]);
    auto tvec = line_point - plane_p0;
    auto qvec = spt::cross(tvec, edges[0]);

    auto t = spt::dot(edges[1], qvec) / spt::dot(edges[0], pvec);
    return line_point + line_dir * t;
}

template <typename Real>
bool does_segment_intersect_triangle(
    const vec<3, Real>& segm_p0, const vec<3, Real>& segm_p1,
    const vec<3, Real>& tr_p0, const vec<3, Real>& tr_p1, const vec<3, Real>& tr_p2) {

    const auto dir = segm_p1 - segm_p0;

    std::array edges{ tr_p1 - tr_p0, tr_p2 - tr_p0 };

    auto pvec = spt::cross(dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    if (std::abs(det) <= std::numeric_limits<Real>::epsilon())
        return false;

    auto inv_det = static_cast<Real>(1) / det;

    auto tvec = segm_p0 - tr_p0;
    auto u = spt::dot(tvec, pvec) * inv_det;
    if (u < static_cast<Real>(0) || u > static_cast<Real>(1))
        return false;

    auto qvec = spt::cross(tvec, edges[0]);
    auto v = spt::dot(dir, qvec) * inv_det;
    if (v < static_cast<Real>(0) || u + v > static_cast<Real>(1))
        return false;

    auto t = spt::dot(edges[1], qvec) * inv_det;
    return t <= static_cast<Real>(1) && t >= static_cast<Real>(0);
}

// todo: return std::optional<vec<3>>
template <typename Real>
bool segment_intersect_plane(
    vec<3, Real>& out_intersect_point,
    const vec<3, Real>& p0, const vec<3, Real>& p1,
    const vec<3, Real>& pl_p0, const vec<3, Real>& pl_p1, const vec<3, Real>& pl_p2) {

    const auto dir = p1 - p0;

    std::array edges{ pl_p1 - pl_p0, pl_p2 - pl_p0 };

    auto pvec = spt::cross(dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    if (std::abs(det) <= std::numeric_limits<Real>::epsilon())
        return false;

    auto tvec = p0 - pl_p0;
    auto qvec = spt::cross(tvec, edges[0]);

    auto t = spt::dot(edges[1], qvec) / det;
    out_intersect_point = p0 + dir * t;
    return t <= static_cast<Real>(1) && t >= static_cast<Real>(0);
}

template <typename Real>
bool does_triangle_intersect_sphere(
    const vec<3, Real>& trngl_p0, const vec<3, Real>& trngl_p1, const vec<3, Real>& trngl_p2,
    const vec<3, Real>& center, Real radius) {

    auto proj = project_on_plane(center, trngl_p0, trngl_p1, trngl_p2);
    if ((proj - center).magnitude2() > radius * radius)
        return false;

    if (is_point_on_triangle(proj, trngl_p0, trngl_p1, trngl_p2, max_sqrs_sum(trngl_p0, trngl_p1, trngl_p2)))
        return true;

    auto closest = closest_triangle_point_to_point_on_plane(proj, trngl_p0, trngl_p1, trngl_p2);
    return (closest - center).magnitude2() <= radius * radius;
}

template <typename Real>
Real sqrs_sum(
    const vec<3, Real>& point,
    const vec<3, Real>& trngl_p0, const vec<3, Real>& trngl_p1, const vec<3, Real>& trngl_p2) {

    std::array sqrs{
        (trngl_p0 - point).magnitude2(),
        (trngl_p1 - point).magnitude2(),
        (trngl_p2 - point).magnitude2()
    };
    return sqrs[0] + sqrs[1] + sqrs[2];
}

template <typename Real>
Real max_sqrs_sum(
    const vec<3, Real>& trngl_p0, const vec<3, Real>& trngl_p1, const vec<3, Real>& trngl_p2) {

    std::array sqrs{
        (trngl_p1 - trngl_p0).magnitude2(),
        (trngl_p2 - trngl_p1).magnitude2(),
        (trngl_p0 - trngl_p2).magnitude2()
    };

    std::array<std::size_t, 2> max_inds;
    if (sqrs[0] < sqrs[1]) {
        max_inds[0] = 1;
        if (sqrs[0] < sqrs[2])
            max_inds[1] = 2;
        else
            max_inds[1] = 0;
    } else {
        max_inds[0] = 0;
        if (sqrs[1] < sqrs[2])
            max_inds[1] = 2;
        else
            max_inds[1] = 1;
    }

    return sqrs[max_inds[0]] + sqrs[max_inds[1]];
}

template <typename Real>
bool is_point_on_triangle(
    const vec<3, Real>& point,
    const vec<3, Real>& trngl_p0, const vec<3, Real>& trngl_p1, const vec<3, Real>& trngl_p2) {

    auto s0 = spt::cross(trngl_p0 - point, trngl_p1 - point).magnitude();
    auto s1 = spt::cross(trngl_p0 - point, trngl_p2 - point).magnitude();
    auto s2 = spt::cross(trngl_p1 - point, trngl_p2 - point).magnitude();
    auto s = spt::cross(trngl_p0 - trngl_p2, trngl_p1 - trngl_p2).magnitude();

    auto expr = s - s0 - s1 - s2;
    return std::abs(expr) <= std::numeric_limits<Real>::epsilon() * (s + s0 + s1 + s2);
}

template <typename Real>
bool is_point_on_triangle(
    const vec<3, Real>& point,
    const vec<3, Real>& trngl_p0, const vec<3, Real>& trngl_p1, const vec<3, Real>& trngl_p2,
    Real max_sqrs_sum) {

    if (sqrs_sum(point, trngl_p0, trngl_p1, trngl_p2) > max_sqrs_sum)
        return false;

    return is_point_on_triangle(point, trngl_p0, trngl_p1, trngl_p2);
}

template <typename Real>
bool is_point_in_tetrahedron(
    const vec<3, Real>& point,
    const vec<3, Real>& tetr_p0, const vec<3, Real>& tetr_p1,
    const vec<3, Real>& tetr_p2, const vec<3, Real>& tetr_p3) {

    auto vert_to_p0 = tetr_p0 - point;
    auto vert_to_p1 = tetr_p1 - point;
    auto vert_to_p2 = tetr_p2 - point;
    auto vert_to_p3 = tetr_p3 - point;

    std::array abs_mixed_prods{
        std::abs(spt::mixed(vert_to_p0, vert_to_p2, vert_to_p3)),
        std::abs(spt::mixed(vert_to_p0, vert_to_p1, vert_to_p2)),
        std::abs(spt::mixed(vert_to_p0, vert_to_p1, vert_to_p3)),
        std::abs(spt::mixed(vert_to_p1, vert_to_p2, vert_to_p3)),
        std::abs(spt::mixed(tetr_p1 - tetr_p0, tetr_p2 - tetr_p0, tetr_p3 - tetr_p0))
    };

    return abs_mixed_prods[0] + abs_mixed_prods[1] + abs_mixed_prods[2] + abs_mixed_prods[3] < abs_mixed_prods[4];
}

template <typename Real>
vec<3, Real> closest_segment_point_to_point(
    const vec<3, Real>& point,
    const vec<3, Real>& segm_p0, const vec<3, Real>& segm_p1) {

    auto proj = project_on_line(point, segm_p0, segm_p1);

    if (weak_in_cuboid(segm_p0, segm_p1, proj)) {
        return proj;
    } else if (std::array sqr_magns{ (segm_p0 - point).magnitude2(), (segm_p1 - point).magnitude2() };
               sqr_magns[0] < sqr_magns[1]) {
        return segm_p0;
    } else {
        return segm_p1;
    }
}

template <typename Real>
vec<3, Real> closest_triangle_point_to_point_on_plane(
    const vec<3, Real>& point,
    const vec<3, Real>& trngl_p0, const vec<3, Real>& trngl_p1, const vec<3, Real>& trngl_p2) {

    std::array closest_points{
        closest_segment_point_to_point(point, trngl_p0, trngl_p1),
        closest_segment_point_to_point(point, trngl_p1, trngl_p2),
        closest_segment_point_to_point(point, trngl_p2, trngl_p0)
    };

    std::array sqrs{
        (closest_points[0] - point).magnitude2(),
        (closest_points[1] - point).magnitude2(),
        (closest_points[2] - point).magnitude2()
    };

    std::size_t min_i;
    if (sqrs[0] < sqrs[1]) {
        if (sqrs[0] < sqrs[2]) {
            min_i = 0;
        } else {
            min_i = 2;
        }
    } else {
        if (sqrs[1] < sqrs[2]) {
            min_i = 1;
        } else {
            min_i = 2;
        }
    }

    return closest_points[min_i];
}

template <typename Real>
Real distance_point_to_line(
    const vec<3, Real>& point,
    const vec<3, Real>& line_p0, const vec<3, Real>& line_p1) {

    return (project_on_line(point, line_p0, line_p1) - point).magnitude();
}

template <typename Real>
Real distance_point_to_segment(
    const vec<3, Real>& point,
    const vec<3, Real>& segm_p0, const vec<3, Real>& segm_p1) {

    auto proj = project_on_line(point, segm_p0, segm_p1);

    if (weak_in_cuboid(segm_p0, segm_p1, proj)) {
        return (proj - point).magnitude();
    } else if (std::array sqr_magns{ (segm_p0 - point).magnitude2(), (segm_p1 - point).magnitude2() };
               sqr_magns[0] < sqr_magns[1]) {
        return std::sqrt(sqr_magns[0]);
    } else {
        return std::sqrt(sqr_magns[1]);
    }
}

template <typename Real>
Real distance_point_to_triangle_on_plane(
    const vec<3, Real>& point,
    const vec<3, Real>& trngl_p0, const vec<3, Real>& trngl_p1, const vec<3, Real>& trngl_p2) {

    std::array closest_points{
        closest_segment_point_to_point(point, trngl_p0, trngl_p1),
        closest_segment_point_to_point(point, trngl_p1, trngl_p2),
        closest_segment_point_to_point(point, trngl_p2, trngl_p0)
    };

    std::array sqrs{
        (closest_points[0] - point).magnitude2(),
        (closest_points[1] - point).magnitude2(),
        (closest_points[2] - point).magnitude2()
    };

    return std::min({ sqrs[0], sqrs[1], sqrs[2] });
}

template <typename Real>
std::pair<vec<3, Real>, vec<3, Real>> lines_closest_points(
    const vec<3, Real>& line0_p0, const vec<3, Real>& line0_p1,
    const vec<3, Real>& line1_p0, const vec<3, Real>& line1_p1) {

    auto u = line0_p1 - line0_p0;
    auto v = line1_p1 - line1_p0;
    auto w = line0_p0 - line1_p0;
    auto a = spt::dot(u, u); // always >= 0
    auto b = spt::dot(u, v);
    auto c = spt::dot(v, v); // always >= 0
    auto d = spt::dot(u, w);
    auto e = spt::dot(v, w);
    auto D = a * c - b * b; // always >= 0

    Real sc, sN, sD = D; // sc = sN / sD, default sD = D >= 0
    Real tc, tN, tD = D; // tc = tN / tD, default tD = D >= 0
    // compute the line parameters of the two closest points
    if (D <= std::numeric_limits<Real>::epsilon()) { // the lines are almost parallel
        sN = static_cast<Real>(0); // force using point P0 on segment S1
        sD = static_cast<Real>(1); // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    } else { // get the closest points on the infinite lines
        sN = b * e - c * d;
        tN = a * e - b * d;
    }

    // finally do the division to get sc and tc
    sc = std::abs(sN) <= std::numeric_limits<Real>::epsilon() ? static_cast<Real>(0) : sN / sD;
    tc = std::abs(tN) <= std::numeric_limits<Real>::epsilon() ? static_cast<Real>(0) : tN / tD;

    return {
        (static_cast<Real>(1) - sc) * line0_p0 + sc * line0_p1,
        (static_cast<Real>(1) - tc) * line1_p0 + tc * line1_p1
    };
}

template <typename Real>
vec<3, Real> lines_closest_point(
    const vec<3, Real>& line0_p0, const vec<3, Real>& line0_p1,
    const vec<3, Real>& line1_p0, const vec<3, Real>& line1_p1) {

    auto closest = lines_closest_points(line0_p0, line0_p1, line1_p0, line1_p1);
    return (std::get<0>(closest) + std::get<1>(closest)) * static_cast<Real>(0.5);
}

template <typename Real>
Real lines_distance(
    const vec<3, Real>& line0_p0, const vec<3, Real>& line0_p1,
    const vec<3, Real>& line1_p0, const vec<3, Real>& line1_p1) {

    auto closest = lines_closest_points(line0_p0, line0_p1, line1_p0, line1_p1);
    auto diff = std::get<0>(closest) - std::get<1>(closest);
    return diff.magnitude();
}

template <typename Real>
std::pair<vec<3, Real>, vec<3, Real>> segments_closest_points(
    const vec<3, Real>& segm0_p0, const vec<3, Real>& segm0_p1,
    const vec<3, Real>& segm1_p0, const vec<3, Real>& segm1_p1) {

    auto u = segm0_p1 - segm0_p0;
    auto v = segm1_p1 - segm1_p0;
    auto w = segm0_p0 - segm1_p0;
    auto a = dot(u, u); // always >= 0
    auto b = dot(u, v);
    auto c = dot(v, v); // always >= 0
    auto d = dot(u, w);
    auto e = dot(v, w);
    auto D = a * c - b * b; // always >= 0
    Real sc, sN, sD = D; // sc = sN / sD, default sD = D >= 0
    Real tc, tN, tD = D; // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D <= std::numeric_limits<Real>::epsilon()) { // the lines are almost parallel
        sN = static_cast<Real>(0); // force using point P0 on segment S1
        sD = static_cast<Real>(1); // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    } else { // get the closest points on the infinite lines
        sN = b * e - c * d;
        tN = a * e - b * d;
        if (sN < static_cast<Real>(0)) { // sc < 0 => the s=0 edge is visible
            sN = static_cast<Real>(0);
            tN = e;
            tD = c;
        } else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < static_cast<Real>(0)) { // tc < 0 => the t=0 edge is visible
        tN = static_cast<Real>(0);
        // recompute sc for this edge
        if (-d < static_cast<Real>(0)) {
            sN = static_cast<Real>(0);
        } else if (-d > a) {
            sN = sD;
        } else {
            sN = -d;
            sD = a;
        }
    } else if (tN > tD) { // tc > 1  => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if (-d + b < static_cast<Real>(0)) {
            sN = static_cast<Real>(0);
        } else if (-d + b > a) {
            sN = sD;
        } else {
            sN = -d + b;
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = std::abs(sN) <= std::numeric_limits<Real>::epsilon() ? static_cast<Real>(0) : sN / sD;
    tc = std::abs(tN) <= std::numeric_limits<Real>::epsilon() ? static_cast<Real>(0) : tN / tD;

    return {
        (static_cast<Real>(1) - sc) * segm0_p0 + sc * segm0_p1,
        (static_cast<Real>(1) - tc) * segm1_p0 + tc * segm1_p1 
    };
}

template <typename Real>
vec<3, Real> segments_closest_point(
    const vec<3, Real>& segm0_p0, const vec<3, Real>& segm0_p1,
    const vec<3, Real>& segm1_p0, const vec<3, Real>& segm1_p1) {

    auto closest = segments_closest_points(segm0_p0, segm0_p1, segm1_p0, segm1_p1);
    return (std::get<0>(closest) + std::get<1>(closest)) * static_cast<Real>(0.5);
}

template <typename Real>
Real segments_distance(
    const vec<3, Real>& segm0_p0, const vec<3, Real>& segm0_p1,
    const vec<3, Real>& segm1_p0, const vec<3, Real>& segm1_p1) {

    auto closest = segments_closest_points(segm0_p0, segm0_p1, segm1_p0, segm1_p1);
    auto diff = std::get<0>(closest) - std::get<1>(closest);
    return diff.magnitude();
}

// CPA - Closest Point of Approach.
template <typename Real>
Real cpa_time(
    const vec<3, Real>& start0, const vec<3, Real>& vel0,
    const vec<3, Real>& start1, const vec<3, Real>& vel1) {

    auto dv = vel0 - vel1;
    auto dv2 = spt::dot(dv, dv);
    if (dv2 <= std::numeric_limits<Real>::epsilon())
        return static_cast<Real>(0);

    auto w0 = start0 - start1;
    return -spt::dot(w0, dv) / dv2;
}

template <typename Real>
Real cpa_distance(
    const vec<3, Real>& start0, const vec<3, Real>& vel0,
    const vec<3, Real>& start1, const vec<3, Real>& vel1) {

    auto time = cpa_time(start0, vel0, start1, vel1);
    auto p0 = start0 + vel0 * time;
    auto p1 = start1 + vel1 * time;
    return (p1 - p0).magnitude();
}

} // namespace spt