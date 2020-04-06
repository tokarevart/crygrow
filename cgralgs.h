// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "vec.h"
#include "mat.h"


namespace cgr {

template <std::size_t Dim>
using pos_t = spt::veci<Dim>;
template <std::size_t Dim>
using upos_t = spt::vecu<Dim>;

template <std::size_t Dim, typename Real = double>
using orientation_t = spt::mat<Dim, Real>; // |each vec| == 1

template <std::size_t Dim>
std::int64_t offset(const pos_t<Dim>& pos, const upos_t<Dim>& dimlens) {
    std::int64_t res = pos.x[0];
    std::int64_t mul = dimlens[0];
    for (std::size_t i = 1; i < Dim; ++i) {
        res += pos.x[i] * mul;
        mul *= dimlens[i];
    }
    return res;
}

template <std::size_t Dim>
upos_t<Dim> upos(std::size_t offset, const upos_t<Dim>& dimlens) {
    upos_t<Dim> res;
    std::size_t t = offset;
    for (std::size_t i = Dim - 1; i > 0; --i) {
        std::size_t prod = dimlens[0];
        for (std::size_t j = 1; j < i; ++j)
            prod *= dimlens[j];
        res[i] = t / prod;
        t -= prod * res[i];
    }
    res[0] = t;
    return res;
}

template <std::size_t Dim>
bool inside(const upos_t<Dim>& pos, const upos_t<Dim>& dimlens) {
    for (std::size_t i = 0; i < Dim; ++i)
        if (pos[i] >= dimlens[i])
            return false;

    return true;
}


template <std::size_t Dim>
using norm_fn = std::function<std::size_t(const pos_t<Dim>&)>;


template <std::size_t Dim, typename ValueType>
std::size_t norm_cryst(
    const pos_t<Dim>& dir,
    const std::vector<spt::vec<Dim, ValueType>>& growdirs) {

    auto cdir = static_cast<spt::vec<Dim, ValueType>>(dir);
    ValueType maxpn = 0;
    for (auto& gd : growdirs) {
        auto pn = std::abs(spt::dot(cdir, gd)) / spt::dot(gd, gd);
        if (pn > maxpn)
            maxpn = pn;
    }
    
    if constexpr (std::is_integral_v<ValueType>)
        return maxpn;
    else
        return maxpn + std::numeric_limits<ValueType>::epsilon();
}

template <std::size_t Dim, typename ValueType>
norm_fn<Dim> make_norm_cryst_fn(std::vector<spt::vec<Dim, ValueType>> growdirs) {
    std::vector<spt::vec<Dim, ValueType>> es_inv_dots;
    es_inv_dots.reserve(growdirs.size());
    for (auto& e : growdirs)
        es_inv_dots.push_back(static_cast<ValueType>(1) / spt::dot(e, e));

    return [std::move(growdirs), std::move(es_inv_dots)](const pos_t<Dim>& dir) -> std::size_t {
        auto cdir = static_cast<spt::vec<Dim, ValueType>>(dir);
        ValueType maxpn = 0;
        for (std::size_t i = 0; i < growdirs.size(); ++i) {
            auto pn = std::abs(spt::dot(cdir, growdirs[i])) * es_inv_dots[i];
            if (pn > maxpn)
                maxpn = pn;
        }
        if constexpr (std::is_integral_v<ValueType>)
            return maxpn;
        else
            return maxpn + std::numeric_limits<ValueType>::epsilon();
    };
}

template <std::size_t Dim>
std::size_t norm_taxicab(const pos_t<Dim>& dir) {
    std::int64_t sum_abs = 0;
    for (auto e : dir.x)
        sum_abs += std::abs(e);
    return static_cast<std::size_t>(sum_abs);
}

template <std::size_t Dim>
std::size_t norm_chebyshev(const pos_t<Dim>& dir) {
    std::size_t max_abs = 0;
    for (auto e : dir.x) {
        auto abse = static_cast<std::size_t>(std::abs(e));
        if (abse > max_abs)
            max_abs = abse;
    }
    return max_abs;
}

template <std::size_t Dim>
std::size_t norm_euclid(const pos_t<Dim>& dir) {
    return std::sqrt(dir.magnitude2()) + std::numeric_limits<double>::epsilon();
}

template <std::size_t Dim>
std::size_t norm2_euclid(const pos_t<Dim>& dir) {
    return dir.magnitude2();
}

} // namespace cgr
