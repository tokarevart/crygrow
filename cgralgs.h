// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "vec.h"


namespace cgr {

template <std::size_t Dim>
std::int64_t offset(const spt::veci<Dim>& pos, const spt::vecu<Dim>& dim_lens) {
    std::int64_t res = pos.x[0];
    std::int64_t mul = dim_lens[0];
    for (std::size_t i = 1; i < Dim; ++i) {
        res += pos.x[i] * mul;
        mul *= dim_lens[i];
    }
    return res;
}

template <std::size_t Dim>
spt::vecu<Dim> upos(std::size_t offset, const spt::vecu<Dim>& dim_lens) {
    spt::vecu<Dim> res;
    std::size_t t = offset;
    for (std::size_t i = Dim - 1; i > 0; --i) {
        std::size_t prod = dim_lens[0];
        for (std::size_t j = 1; j < i; ++j)
            prod *= dim_lens[j];
        res[i] = t / prod;
        t -= prod * res[i];
    }
    res[0] = t;
    return res;
}

template <std::size_t Dim>
bool inside(const spt::vecu<Dim>& pos, const spt::vecu<Dim>& dim_lens) {
    for (std::size_t i = 0; i < Dim; ++i)
        if (pos[i] >= dim_lens[i])
            return false;

    return true;
}


template <std::size_t Dim>
using norm_fn = std::function<std::size_t(const spt::veci<Dim>&)>;


template <std::size_t Dim, typename ValueType>
std::size_t norm_cryst(
    const spt::veci<Dim>& dir, 
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

    return [std::move(growdirs), std::move(es_inv_dots)](const spt::veci<Dim>& dir) -> std::size_t {
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
std::size_t norm_taxicab(const spt::veci<Dim>& dir) {
    std::int64_t sum_abs = 0;
    for (auto e : dir.x)
        sum_abs += std::abs(e);
    return static_cast<std::size_t>(sum_abs);
}

template <std::size_t Dim>
std::size_t norm_chebyshev(const spt::veci<Dim>& dir) {
    std::size_t max_abs = 0;
    for (auto e : dir.x) {
        auto abse = static_cast<std::size_t>(std::abs(e));
        if (abse > max_abs)
            max_abs = abse;
    }
    return max_abs;
}

template <std::size_t Dim>
std::size_t norm_euclid(const spt::veci<Dim>& dir) {
    return std::sqrt(dir.magnitude2()) + std::numeric_limits<double>::epsilon();
}

template <std::size_t Dim>
std::size_t norm2_euclid(const spt::veci<Dim>& dir) {
    return dir.magnitude2();
}

} // namespace cgr
