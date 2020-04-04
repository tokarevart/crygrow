// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <vector>
#include <functional>
#include "vec.h"
#include "cgralgs.h"


namespace cgr {

enum class nbhood_kind {
    von_neumann,
    moore,
    euclid,
    crystallographic
};


template <std::size_t Dim>
bool inside_nbhood(norm_fn<Dim> normfn, const spt::veci<Dim>& pos, std::size_t range) {
    return normfn(pos) <= range;
}


using nbhood_offset = std::vector<std::size_t>;

template <std::size_t Dim, typename Cell>
using get_cell_by_offset_fn = std::function<Cell*(const std::size_t)>;

template <std::size_t Dim, typename Cell>
using try_get_cell_fn = std::function<Cell*(const spt::veci<Dim>&)>;

template <std::size_t Dim>
using inside_fn = std::function<bool(const spt::veci<Dim>&)>;


template <std::size_t Dim>
using nbhood_pos = std::vector<spt::veci<Dim>>;

template <std::size_t Dim>
using nbhood_pos_shift_fn = std::function<nbhood_pos<Dim>(std::size_t)>;

template <std::size_t Dim, typename Real = double>
using orientation_t = spt::mat<Dim, Real>; // |each vec| == 1


template <std::size_t Dim>
nbhood_pos<Dim> make_nbhood_pos_shift(norm_fn<Dim> normfn, std::size_t range) {
    std::vector<spt::veci<Dim>> res;
    std::size_t buf = 2 * range + 1;
    std::int64_t srange = range;
    if constexpr (Dim == 2) {
        res.reserve(buf * buf - 1);
        for (std::int64_t y = -srange; y <= srange; ++y) {
            for (std::int64_t x = -srange; x <= srange; ++x) {
                auto sh = spt::veci<Dim>{ x, y };
                if (!(x == 0 && y == 0) &&
                    inside_nbhood<Dim>(normfn, sh, range))
                    res.push_back(sh);
            }
        }
        res.shrink_to_fit();

    } else if constexpr (Dim == 3) {
        res.reserve(buf * buf * buf - 1);
        for (std::int64_t z = -srange; z <= srange; ++z) {
            for (std::int64_t y = -srange; y <= srange; ++y) {
                for (std::int64_t x = -srange; x <= srange; ++x) {
                    auto sh = spt::veci<Dim>{ x, y, z };
                    if (!(x == 0 && y == 0 && z == 0) &&
                        inside_nbhood<Dim>(normfn, sh, range))
                        res.push_back(sh);
                }
            }
        }
        res.shrink_to_fit();

    } else {
        static_assert(false);
    }

    return res;
}


template <std::size_t Dim>
nbhood_pos_shift_fn<Dim> make_nbhood_pos_shift_fn(norm_fn<Dim> normfn) {
    return [normfn](std::size_t range) -> nbhood_pos<Dim> {
        return make_nbhood_pos_shift<Dim>(normfn, range);
    };
}


template <std::size_t Dim>
nbhood_offset nbhood_pos_to_offset(const nbhood_pos<Dim>& nbhpos, const spt::vecu<Dim>& dim_lens) {
    nbhood_offset res;
    res.reserve(nbhpos.size());
    for (auto pos : nbhpos)
        res.push_back(cgr::offset(pos, dim_lens));
    return res;
}


template <std::size_t Dim>
nbhood_pos<Dim> apply_nbhood_pos_shift(
    const spt::veci<Dim>& center, const nbhood_pos<Dim>& nbhshift,
    std::optional<inside_fn<Dim>> infn) {

    nbhood_pos<Dim> res;
    res.reserve(nbhshift.size());
    for (auto shift : nbhshift) {
        auto new_pos = center + shift;
        if (!infn)
            res.push_back(new_pos);
        else 
            if (infn.value()(new_pos))
                res.push_back(new_pos);
    }
        
    return res;
}


template <std::size_t Dim>
nbhood_pos<Dim> make_nbhood_pos(
    nbhood_pos_shift_fn<Dim> shfn, const spt::veci<Dim>& center, std::size_t range,
    std::optional<inside_fn<Dim>> infn) {

    return apply_nbhood_pos_shift<Dim>(center, shfn(range), infn);
}


template <std::size_t Dim>
nbhood_pos<Dim> make_nbhood_pos(
    norm_fn<Dim> normfn, const spt::veci<Dim>& center, std::size_t range,
    std::optional<inside_fn<Dim>> infn) {

    return apply_nbhood_pos_shift<Dim>(center, make_nbhood_pos_shift<Dim>(normfn, range), infn);
}


template <std::size_t Dim>
nbhood_offset make_nbhood_offset(
    nbhood_pos_shift_fn<Dim> shfn, std::size_t center, const spt::vecu<Dim>& dim_lens,
    std::size_t range, std::optional<inside_fn<Dim>> infn) {

    spt::veci<Dim> center_pos = cgr::upos(center, dim_lens);
    auto nbhpos = make_nbhood_pos<Dim>(shfn, center_pos, range, infn);
    return nbhood_pos_to_offset<Dim>(nbhpos, dim_lens);
}


template <std::size_t Dim>
nbhood_offset make_nbhood_offset(
    norm_fn<Dim> normfn, std::size_t center, const spt::vecu<Dim>& dim_lens,
    std::size_t range, std::optional<inside_fn<Dim>> infn) {

    spt::veci<Dim> center_pos = cgr::upos(center, dim_lens);
    auto nbhpos = make_nbhood_pos<Dim>(normfn, center_pos, range, infn);
    return nbhood_pos_to_offset<Dim>(nbhpos, dim_lens);
}


template <std::size_t Dim, typename Cell>
using nbhood = std::vector<Cell*>;


template <std::size_t Dim, typename Cell>
nbhood<Dim, Cell> make_nbhood(
    norm_fn<Dim> normfn, const spt::veci<Dim>& center, std::size_t range,
    try_get_cell_fn<Dim, Cell> trygetcell) {

    nbhood<Dim, Cell> res;
    auto nbhpos = trygetcell ? 
        make_nbhood_pos<Dim>(normfn, center, range,
            [trygetcell](const spt::veci<Dim>& pos) -> bool { 
                return trygetcell.value()(pos); }) 
        : make_nbhood_pos<Dim>(normfn, center, range);

    for (auto& pos : nbhpos) {
        auto pcell = trygetcell(pos);
        if (pcell)
            res.push_back(pcell);
    }

    return res;
}


template <std::size_t Dim, typename Cell>
nbhood<Dim, Cell> make_nbhood(
    const nbhood_offset& nbhoff,
    get_cell_by_offset_fn<Dim, Cell> getcell) {

    nbhood<Dim, Cell> res;
    res.reserve(nbhoff.size());
    for (auto& off : nbhoff)
        res.push_back(getcell(off));
    return res;
}


template <std::size_t Dim, typename Cell>
nbhood<Dim, Cell> make_nbhood(
    const nbhood_pos<Dim>& nbhpos,
    try_get_cell_fn<Dim, Cell> trygetcell) {
    nbhood<Dim, Cell> res;
    for (auto& pos : nbhpos) {
        auto pcell = trygetcell(pos);
        if (pcell)
            res.push_back(pcell);
    }

    return res;
}

} // namespace cgr
