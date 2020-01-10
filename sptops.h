// Copyright © 2018-2020 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cmath>
#include <algorithm>
#include "mat.h"

// todo: try vectorization (SIMD) and make benchmark
namespace spt {

template <std::size_t Dim, typename ValueType>
mat<Dim, ValueType> dot(const mat<Dim, ValueType>& mat0, const mat<Dim, ValueType>& mat1) {
    mat<Dim, ValueType> res;
    for (std::size_t i = 0; i < Dim; ++i) {
        for (std::size_t j = 0; j < Dim; ++j) {
            ValueType buf = 0.0;
            for (std::size_t k = 0; k < Dim; ++k)
                buf += mat0[i][k] * mat1[k][j];
            res[i][j] = buf;
        }
    }
    return res;
}

template <std::size_t Dim, typename ValueType>
vec<Dim, ValueType> dot(const mat<Dim, ValueType>& matr, const vec<Dim, ValueType>& vect) {
    vec<Dim, ValueType> res;
    for (std::size_t i = 0; i < Dim; ++i)
        res[i] = dot(matr[i], vect);
    return res;
}

template <std::size_t Dim, typename ValueType>
ValueType dot(const vec<Dim, ValueType>& vect0, const vec<Dim, ValueType>& vect1) {
    ValueType res = 0;
    for (std::size_t i = 0; i < Dim; ++i)
        res += vect0[i] * vect1[i];
    return res;
}

template <typename ValueType>
vec<3, ValueType> cross(const vec<3, ValueType>& vect0, const vec<3, ValueType>& vect1) {
    return {
        vect0.x[1] * vect1.x[2] - vect0.x[2] * vect1.x[1],
        vect0.x[2] * vect1.x[0] - vect0.x[0] * vect1.x[2],
        vect0.x[0] * vect1.x[1] - vect0.x[1] * vect1.x[0] };
}

template <typename ValueType>
ValueType cross(const vec<2, ValueType>& vect0, const vec<2, ValueType>& vect1) {
    return vect0[0] * vect1[1] - vect0[1] * vect1[0];
}

template <typename ValueType>
ValueType mixed(const vec<3, ValueType>& vect0, 
                const vec<3, ValueType>& vect1, 
                const vec<3, ValueType>& vec2) {
    return dot(cross(vect0, vect1), vec2);
}

template <std::size_t Dim, typename Real>
Real cos(const vec<Dim, Real>& vect0, const vec<Dim, Real>& vect1) {
    return dot(vect0, vect1) / std::sqrt(vect0.magnitude2() * vect1.magnitude2());
}

template <std::size_t Dim, typename ValueType>
void sort_elementwise(vec<Dim, ValueType>& vect0, vec<Dim, ValueType>& vect1) {
    for (std::size_t i = 0; i < Dim; ++i)
        if (vect0[i] > vect1[i])
            std::swap(vect0[i], vect1[i]);
}

} // namespace spt