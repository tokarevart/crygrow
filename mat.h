// Copyright © 2018-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "vec.h"

namespace spt {

using default_value_type = double;

template <std::size_t N, typename ValueType = default_value_type>
struct mat;

template <std::size_t N> using mati = mat<N, std::int64_t>;
template <std::size_t N> using matf = mat<N, float>;
template <std::size_t N> using matd = mat<N, double>;
template <typename ValueType> 
using mat2 = mat<2, ValueType>;
using mat2i = mati<2>;
using mat2f = matf<2>;
using mat2d = matd<2>;
template <typename ValueType> 
using mat3 = mat<3, ValueType>;
using mat3i = mati<3>;
using mat3f = matf<3>;
using mat3d = matd<3>;
template <typename ValueType> 
using mat4 = mat<4, ValueType>;
using mat4i = mati<4>;
using mat4f = matf<4>;
using mat4d = matd<4>;

template <typename ValueType>
struct mat<3, ValueType> {
    static constexpr std::size_t n = 3;
    using value_type = ValueType;
    using line_type = spt::vec<3, value_type>;

    std::array<line_type, 3> x;

    static mat identity() {
        return { 
            line_type({ 
            static_cast<value_type>(1), 
            static_cast<value_type>(0), 
            static_cast<value_type>(0) }),
            line_type({ 
            static_cast<value_type>(0), 
            static_cast<value_type>(1), 
            static_cast<value_type>(0) }),
            line_type({ 
            static_cast<value_type>(0), 
            static_cast<value_type>(0), 
            static_cast<value_type>(1) }) };
    }

    mat& transpose() const {
        std::swap(x[0][1], x[1][0]);
        std::swap(x[0][2], x[2][0]);
        std::swap(x[1][2], x[2][1]);
        return *this;
    }
    mat transposed() const {
        return {
            x[0][0], x[1][0], x[2][0],
            x[0][1], x[1][1], x[2][1],
            x[0][2], x[1][2], x[2][2] };
    }
    mat inversed() const {
        auto a00 = x[0][0], a01 = x[0][1], a02 = x[0][2];
        auto a10 = x[1][0], a11 = x[1][1], a12 = x[1][2];
        auto a20 = x[2][0], a21 = x[2][1], a22 = x[2][2];

        auto b01 = a22 * a11 - a12 * a21;
        auto b11 = -a22 * a10 + a12 * a20;
        auto b21 = a21 * a10 - a11 * a20;

        auto det = a00 * b01 + a01 * b11 + a02 * b21;

        return mat{ b01, (-a22 * a01 + a02 * a21), (a12 * a01 - a02 * a11),
                    b11, (a22 * a00 - a02 * a20), (-a12 * a00 + a02 * a10),
                    b21, (-a21 * a00 + a01 * a20), (a11 * a00 - a01 * a10) } / det;
    }

    mat& operator=(const mat& right) {
        x = right.x;
        return *this;
    }
    mat operator-() const {
        return { -x[0], -x[1], -x[2] };
    }
    mat operator+(const mat& right) const {
        return {
            x[0] + right.x[0],
            x[1] + right.x[1],
            x[2] + right.x[2] };
    }
    mat operator-(const mat& right) const {
        return {
            x[0] - right.x[0],
            x[1] - right.x[1],
            x[2] - right.x[2] };
    }
    mat operator*(value_type scalar) const {
        return {
            x[0] * scalar,
            x[1] * scalar,
            x[2] * scalar };
    }
    mat operator/(value_type scalar) const {
        return {
            x[0] / scalar,
            x[1] / scalar,
            x[2] / scalar };
    }
    mat& operator+=(const mat& right) {
        x[0] += right.x[0];
        x[1] += right.x[1];
        x[2] += right.x[2];
        return *this;
    }
    mat& operator-=(const mat& right) {
        x[0] -= right.x[0];
        x[1] -= right.x[1];
        x[2] -= right.x[2];
        return *this;
    }
    mat& operator*=(value_type scalar) {
        x[0] *= scalar;
        x[1] *= scalar;
        x[2] *= scalar;
        return *this;
    }
    mat& operator/=(value_type scalar) {
        x[0] /= scalar;
        x[1] /= scalar;
        x[2] /= scalar;
        return *this;
    }
    line_type& operator[](std::size_t i) {
        return x[i];
    }
    const line_type& operator[](std::size_t i) const {
        return x[i];
    }

    mat() {}
    mat(const mat& other) {
        x = other.x;
    }
    mat(const std::array<line_type, 3>& x)
        : x{ x } {}
    mat(const std::array<value_type, 9>& x) {
        this->x[0] = line_type({ x[0], x[1], x[2] });
        this->x[1] = line_type({ x[3], x[4], x[5] });
        this->x[2] = line_type({ x[6], x[7], x[8] });
    }
    mat(const line_type& x0, const line_type& x1, const line_type& x2) {
        x[0] = x0;
        x[1] = x1;
        x[2] = x2;
    }
    mat(value_type v0, value_type v1, value_type v2,
        value_type v3, value_type v4, value_type v5,
        value_type v6, value_type v7, value_type v8)
        : x{ line_type{ v0, v1, v2 }, line_type{ v3, v4, v5 }, line_type{ v6, v7, v8 } } {}
};

template <typename ValueType>
struct mat<2, ValueType> {
    static constexpr std::size_t n = 2;
    using value_type = ValueType;
    using line_type = spt::vec<2, value_type>;

    std::array<line_type, 2> x;

    static mat identity() {
        return {
            line_type({ static_cast<value_type>(1), static_cast<value_type>(0) }),
            line_type({ static_cast<value_type>(0), static_cast<value_type>(1) }) };
    }

    mat& transpose() const {
        std::swap(x[0][1], x[1][0]);
        return *this;
    }
    mat transposed() const {
        return {
            x[0][0], x[1][0],
            x[0][1], x[1][1] };
    }
    mat inversed() const {
        value_type det = x[0][0] * x[1][1] - x[0][1] * x[1][0];
        return mat{ x[1][1], -x[0][1],
                   -x[1][0],  x[0][0] } / det;
    }

    mat& operator=(const mat& right) {
        x = right.x;
        return *this;
    }
    mat operator-() const {
        return { -x[0], -x[1] };
    }
    mat operator+(const mat& right) const {
        return {
            x[0] + right.x[0],
            x[1] + right.x[1] };
    }
    mat operator-(const mat& right) const {
        return {
            x[0] - right.x[0],
            x[1] - right.x[1] };
    }
    mat operator*(value_type scalar) const {
        return {
            x[0] * scalar,
            x[1] * scalar };
    }
    mat operator/(value_type scalar) const {
        return {
            x[0] / scalar,
            x[1] / scalar };
    }
    mat& operator+=(const mat& right) {
        x[0] += right.x[0];
        x[1] += right.x[1];
        return *this;
    }
    mat& operator-=(const mat& right) {
        x[0] -= right.x[0];
        x[1] -= right.x[1];
        return *this;
    }
    mat& operator*=(value_type scalar) {
        x[0] *= scalar;
        x[1] *= scalar;
        return *this;
    }
    mat& operator/=(value_type scalar) {
        x[0] /= scalar;
        x[1] /= scalar;
        return *this;
    }
    line_type& operator[](std::size_t i) {
        return x[i];
    }
    const line_type& operator[](std::size_t i) const {
        return x[i];
    }

    mat() {}
    mat(const mat& other) {
        x = other.x;
    }
    mat(const std::array<line_type, 2>& x)
        : x{x} {}
    mat(const std::array<value_type, 4>& x) {
        this->x[0] = line_type({ x[0], x[1] });
        this->x[1] = line_type({ x[2], x[3] });
    }
    mat(const line_type& x0, const line_type& x1) {
        x[0] = x0;
        x[1] = x1;
    }
    mat(value_type v0, value_type v1,
        value_type v2, value_type v3)
        : x{ line_type({ v0, v1 }), line_type({ v2, v3 }) } {}
};

}  // namespace spt