#pragma once
#include "vec.h"


namespace cgr {

template <std::size_t Dim>
std::size_t offset(const spt::vecu<Dim>& pos, const spt::vecu<Dim>& dim_lens) {
    std::size_t res = pos.x[0];
    std::size_t mul = dim_lens[0];
    for (std::size_t i = 1; i < Dim; ++i) {
        res += pos.x[i] * mul;
        mul *= dim_lens[i];
    }
    return res;
}

template <std::size_t Dim>
spt::vecu<Dim> upos(std::size_t offset, const spt::vecu<Dim>& dim_lens) {
    spt::vecu<Dim> res;
    if constexpr (Dim == 2) {
        auto x = dim_lens[0];
        res[1] = offset / x;
        res[0] = offset - x * res[1];

    } else if constexpr (Dim == 3) {
        auto x = dim_lens[0];
        auto y = dim_lens[1];
        auto xy = x * y;
        res[2] = offset / xy;
        auto t = offset - xy * res[2];
        res[1] = t / x;
        res[0] = t - x * res[1];

    } else {
        static_assert(false);
    }

    return res;
}

template <std::size_t Dim>
bool inside(const spt::vecu<Dim>& pos, const spt::vecu<Dim>& dim_lens) {
    for (std::size_t i = 0; i < Dim; ++i)
        if (pos[i] >= dim_lens[i])
            return false;

    return true;
}

} // namespace cgr
