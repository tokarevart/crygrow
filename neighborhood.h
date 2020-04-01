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
    euclid
};

template <std::size_t Dim>
bool inside_nbhood(nbhood_kind kind, const spt::veci<Dim>& pos, std::size_t range) {
    switch (kind) {
    case nbhood_kind::von_neumann:
        return cgr::norm_taxicab(pos) <= range;

    case nbhood_kind::moore:
        for (auto e : pos.x)
            if (e > static_cast<std::int64_t>(range))
                return false;
        return true;

    case nbhood_kind::euclid:
        return cgr::norm2_euclid(pos) <= range * range;
    }

    std::terminate();
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
nbhood_pos<Dim> make_nbhood_pos(
    nbhood_kind kind, const spt::veci<Dim>& center, std::size_t range,
    std::optional<inside_fn<Dim>> infn) {

    std::int64_t srange = range;
    std::vector<spt::veci<Dim>> res;
    if constexpr (Dim == 2) {
        res.reserve(4 * range * (range + 1)); // (2 * range + 1)^2 - 1
        for (std::int64_t y = -srange; y <= srange; ++y) {
            for (std::int64_t x = -srange; x <= srange; ++x) {
                auto rel = spt::veci<Dim>{ x, y };
                if (!(x == 0 && y == 0) &&
                    inside_nbhood<Dim>(kind, rel, range)) {
                    auto new_pos = center + rel;
                    if (!infn)
                        res.push_back(new_pos);
                    else
                        if (infn.value()(new_pos))
                            res.push_back(new_pos);
                }
            }
        }
                
        res.shrink_to_fit();

    } else if constexpr (Dim == 3) {
        std::size_t buf = 2 * range + 1;
        res.reserve(buf * buf * buf - 1);
        for (std::int64_t z = -srange; z <= srange; ++z) {
            for (std::int64_t y = -srange; y <= srange; ++y) {
                for (std::int64_t x = -srange; x <= srange; ++x) {
                    auto rel = spt::veci<Dim>{ x, y, z };
                    if (!(x == 0 && y == 0 && z == 0) &&
                        inside_nbhood<Dim>(kind, rel, range)) {
                        auto new_pos = center + rel;
                        if (!infn)
                            res.push_back(new_pos);
                        else
                            if (infn.value()(new_pos))
                                res.push_back(new_pos);
                    }
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
nbhood_pos<Dim> make_nbhood_pos_shift(nbhood_kind kind, std::size_t range) {
    return make_nbhood_pos<Dim>(kind, spt::veci<Dim>::filled_with(0), range, std::nullopt);
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

    nbhood_offset res;
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
nbhood_offset make_nbhood_offset(
    nbhood_kind kind, std::size_t center, const spt::vecu<Dim>& dim_lens, std::size_t range,
    std::optional<inside_fn<Dim>> infn) {

    spt::veci<Dim> center_pos = cgr::upos(center, dim_lens);
    auto nbhpos = make_nbhood_pos<Dim>(kind, center_pos, range, infn);
    return nbhood_pos_to_offset<Dim>(nbhpos, dim_lens);
}


template <std::size_t Dim, typename Cell>
using nbhood = std::vector<Cell*>;


template <std::size_t Dim, typename Cell>
nbhood<Dim, Cell> make_nbhood(
    nbhood_kind kind, const spt::veci<Dim>& center, std::size_t range,
    try_get_cell_fn<Dim, Cell> trygetcell) {

    nbhood<Dim, Cell> res;
    auto nbhpos = trygetcell ? 
        make_nbhood_pos<Dim>(kind, center, range,
            [trygetcell](const spt::veci<Dim>& pos) -> bool { 
                return trygetcell.value()(pos); }) 
        : make_nbhood_pos<Dim>(kind, center, range);

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
