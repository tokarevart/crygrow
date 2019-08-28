// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cmath>
#include <cstdint>
#include <array>

namespace spt {

using default_value_type = double;

template <std::size_t Dim, typename ValueType = default_value_type>
struct vec;

template <std::size_t Dim> using veci = vec<Dim, std::int64_t>;
template <std::size_t Dim> using vecf = vec<Dim, float>;
template <std::size_t Dim> using vecd = vec<Dim, double>;
template <typename ValueType> 
using vec2 = vec<2, ValueType>;
using vec2i = vec2<std::int64_t>;
using vec2f = vec2<float>;
using vec2d = vec2<double>;
template <typename ValueType> 
using vec3 = vec<3, ValueType>;
using vec3i = vec3<std::int64_t>;
using vec3f = vec3<float>;
using vec3d = vec3<double>;
template <typename ValueType> 
using vec4 = vec<4, ValueType>;
using vec4i = vec4<std::int64_t>;
using vec4f = vec4<float>;
using vec4d = vec4<double>;


template <std::size_t Dim, typename ValueType>
struct vec {
    static constexpr std::size_t dim = Dim;
    using value_type = ValueType;

    std::array<ValueType, Dim> x;

    ValueType magnitude() const {
        return std::sqrt(magnitude2());
    }
    ValueType magnitude2() const {
        return dot(*this, *this);
    }

    vec& normalize() {
        auto inv_magn = static_cast<ValueType>(1) / magnitude();
        //x[0] *= inv_magn;
        //x[1] *= inv_magn;
        //x[2] *= inv_magn;
        return *this;
    }

    vec& operator=(const vec& right) {
        x = right.x;
        return *this;
    }
    vec operator-() const {
        //return { -x[0], -x[1], -x[2] };
    }
    vec operator+(const vec& right) const {
        //return {
        //    x[0] + right.x[0],
        //    x[1] + right.x[1],
        //    x[2] + right.x[2] };
    }
    vec operator-(const vec& right) const {
        //return {
        //    x[0] - right.x[0],
        //    x[1] - right.x[1],
        //    x[2] - right.x[2] };
    }
    vec operator*(ValueType scalar) const {
        //return {
        //    x[0] * scalar,
        //    x[1] * scalar,
        //    x[2] * scalar };
    }
    vec operator/(ValueType scalar) const {
        //return {
        //    x[0] / scalar,
        //    x[1] / scalar,
        //    x[2] / scalar };
    }
    vec& operator+=(const vec& right) {
        //x[0] += right.x[0];
        //x[1] += right.x[1];
        //x[2] += right.x[2];
        return *this;
    }
    vec& operator-=(const vec& right) {
        //x[0] -= right.x[0];
        //x[1] -= right.x[1];
        //x[2] -= right.x[2];
        return *this;
    }
    vec& operator*=(ValueType scalar) {
        //x[0] *= scalar;
        //x[1] *= scalar;
        //x[2] *= scalar;
        return *this;
    }
    vec& operator/=(ValueType scalar) {
        //x[0] /= scalar;
        //x[1] /= scalar;
        //x[2] /= scalar;
        return *this;
    }
    ValueType& operator[](std::size_t i) {
        return x[i];
    }
    const ValueType& operator[](std::size_t i) const {
        return x[i];
    }
    template <typename NewValueType>
    operator vec<3, NewValueType>() const {
        //return {
        //    static_cast<NewValueType>(x[0]),
        //    static_cast<NewValueType>(x[1]),
        //    static_cast<NewValueType>(x[2])
        //};
    }

    vec() {
        std::fill(x.begin(), x.end(), static_cast<ValueType>(0));
    }
    vec(const vec& other) {
        x = other.x;
    }
    vec(const std::array<ValueType, Dim> & x)
        : x{ x } {}
    template <typename... Values>
    vec(Values&&... xs)
        : x{ std::forward<ValueType>(static_cast<ValueType>(xs))... } {}
};


template<std::size_t Dim, typename ValueType>
vec(const std::array<ValueType, Dim>& x) -> vec<Dim, ValueType>;

template<typename ValueType, typename... TaleValues>
vec(ValueType x0, TaleValues... xt) -> vec<1 + sizeof...(TaleValues), ValueType>;


template<std::size_t Dim, typename ValueType>
struct std::hash<spt::vec<Dim, ValueType>> {
    std::size_t operator() (const spt::vec<Dim, ValueType>& key) const {
        std::hash<ValueType> hasher;
        std::size_t h = 0;
        for (auto e : key.x)
            h ^= hasher(e) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

} // namespace spt


template <std::size_t Dim, typename ValueType>
spt::vec<Dim, ValueType> operator*(ValueType scalar, const spt::vec<Dim, ValueType>& v) {
    return v * scalar;
}
