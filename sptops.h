// Copyright © 2018-2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cmath>
#include <algorithm>
#include "mat.h"

// TODO: try vectorization (SIMD) and make benchmark
namespace spt {

template <std::size_t Dim, typename ValueType>
mat<Dim, ValueType> dot(const mat<Dim, ValueType>& mat0, const mat<Dim, ValueType>& mat1) {
    //std::array mat1_cols{
    //    vec<3, Real>{ mat1[0][0], mat1[1][0], mat1[2][0] },
    //    vec<3, Real>{ mat1[0][1], mat1[1][1], mat1[2][1] },
    //    vec<3, Real>{ mat1[0][2], mat1[1][2], mat1[2][2] } };
    //return {
    //    { dot(mat0[0], mat1_cols[0]), dot(mat0[0], mat1_cols[1]), dot(mat0[0], mat1_cols[2]) },
    //    { dot(mat0[1], mat1_cols[0]), dot(mat0[1], mat1_cols[1]), dot(mat0[1], mat1_cols[2]) },
    //    { dot(mat0[2], mat1_cols[0]), dot(mat0[2], mat1_cols[1]), dot(mat0[2], mat1_cols[2]) } };
}

template <std::size_t Dim, typename ValueType>
vec<Dim, ValueType> dot(const mat<Dim, ValueType>& matr, const vec<Dim, ValueType>& vect) {
    //return { dot(matr[0], vect), dot(matr[1], vect), dot(matr[2], vect) };
}

template <std::size_t Dim, typename ValueType>
ValueType dot(const vec<Dim, ValueType>& vec0, const vec<Dim, ValueType>& vec1) {
    //return vec0[0] * vec1[0] + vec0[1] * vec1[1] + vec0[2] * vec1[2];
}

template <typename ValueType>
vec<3, ValueType> cross(const vec<3, ValueType>& vec0, const vec<3, ValueType>& vec1) {
    return {
        vec0.x[1] * vec1.x[2] - vec0.x[2] * vec1.x[1],
        vec0.x[2] * vec1.x[0] - vec0.x[0] * vec1.x[2],
        vec0.x[0] * vec1.x[1] - vec0.x[1] * vec1.x[0] };
}

template <typename ValueType>
ValueType cross(const vec<2, ValueType>& vec0, const vec<2, ValueType>& vec1) {
    return vec0[0] * vec1[1] - vec0[1] * vec1[0];
}

template <typename ValueType>
ValueType mixed(const vec<3, ValueType>& vec0, 
                const vec<3, ValueType>& vec1, 
                const vec<3, ValueType>& vec2) {
    return dot(cross(vec0, vec1), vec2);
}

template <std::size_t Dim, typename Real>
Real cos(const vec<Dim, Real>& vec0, const vec<Dim, Real>& vec1) {
    return dot(vec0, vec1) / std::sqrt(vec0.magnitude2() * vec1.magnitude2());
}

} // namespace spt