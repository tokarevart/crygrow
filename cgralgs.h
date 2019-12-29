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
std::size_t magnitude_von_neumann(const spt::veci<Dim>& pos) {
    std::int64_t sum_abs = 0;
    for (auto e : pos.x)
        sum_abs += std::abs(e);
    return static_cast<std::size_t>(sum_abs);
}

template <std::size_t Dim>
std::size_t magnitude_moore(const spt::veci<Dim>& pos) {
    std::size_t max_abs = 0;
    for (auto e : pos.x) {
        auto abse = static_cast<std::size_t>(std::abs(e));
        if (abse > max_abs)
            max_abs = abse;
    }
    return max_abs;
}

template <std::size_t Dim, typename Real = double>
Real magnitude_euclid(const spt::veci<Dim>& pos) {
    return std::sqrt(pos.magnitude2());
}

template <std::size_t Dim>
std::size_t magnitude2_euclid(const spt::veci<Dim>& pos) {
    return pos.magnitude2();
}

} // namespace cgr
