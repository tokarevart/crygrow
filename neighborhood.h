// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <vector>
#include <functional>
#include "vec.h"
#include "cgralgs.h"


namespace cgr::nbh {

enum class nbhood_kind {
    von_neumann,
    moore,
    euclid,
    crystallographic
};


template <std::size_t Dim>
bool inside_nbhood(norm_fn<Dim> normfn, const pos_t<Dim>& pos, std::size_t range) {
    return normfn(pos) <= range;
}


template <std::size_t Dim>
using shifts_fn = std::function<std::vector<pos_t<Dim>>(std::size_t)>;

using offsets_t = std::vector<std::size_t>;

template <std::size_t Dim, typename Cell>
using cell_by_offset_fn = std::function<Cell*(const std::size_t)>;

template <std::size_t Dim, typename Cell>
using try_cell_fn = std::function<Cell*(const pos_t<Dim>&)>;

template <std::size_t Dim>
using inside_fn = std::function<bool(const pos_t<Dim>&)>;


template <std::size_t Dim>
std::vector<pos_t<Dim>> make_shifts(norm_fn<Dim> normfn, std::size_t range, std::size_t bbox_range = 0) {
    if (bbox_range == 0)
        bbox_range = range;

    std::vector<pos_t<Dim>> res;
    std::size_t buf = 2 * bbox_range + 1;
    std::int64_t sbbox_range = bbox_range;
    if constexpr (Dim == 2) {
        res.reserve(buf * buf - 1);
        for (std::int64_t y = -sbbox_range; y <= sbbox_range; ++y) {
            for (std::int64_t x = -sbbox_range; x <= sbbox_range; ++x) {
                auto sh = pos_t<Dim>({ x, y });
                if (!(x == 0 && y == 0) &&
                    inside_nbhood<Dim>(normfn, sh, range))
                    res.push_back(sh);
            }
        }
        res.shrink_to_fit();

    } else if constexpr (Dim == 3) {
        res.reserve(buf * buf * buf - 1);
        for (std::int64_t z = -sbbox_range; z <= sbbox_range; ++z) {
            for (std::int64_t y = -sbbox_range; y <= sbbox_range; ++y) {
                for (std::int64_t x = -sbbox_range; x <= sbbox_range; ++x) {
                    auto sh = pos_t<Dim>({ x, y, z });
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
shifts_fn<Dim> make_shifts_fn(norm_fn<Dim> normfn) {
    return [normfn](std::size_t range) -> std::vector<pos_t<Dim>> {
        return make_shifts<Dim>(normfn, range);
    };
}


template <std::size_t Dim>
offsets_t poses_to_offsets(const std::vector<pos_t<Dim>>& poses, const spt::vecu<Dim>& dimlens) {
    offsets_t res;
    res.reserve(poses.size());
    for (auto pos : poses)
        res.push_back(cgr::offset(pos, dimlens));
    return res;
}

template <std::size_t Dim, typename InsideFn>
std::vector<pos_t<Dim>> apply_shifts(
    const pos_t<Dim>& pos, const std::vector<pos_t<Dim>>& shifts,
    InsideFn infn) {

    std::vector<pos_t<Dim>> res;
    res.reserve(shifts.size());
    for (auto shift : shifts) {
        auto new_pos = pos + shift;
        if (infn(new_pos))
            res.push_back(new_pos);
    }
        
    return res;
}

template <std::size_t Dim>
std::vector<pos_t<Dim>> apply_shifts(
    const pos_t<Dim>& pos, const std::vector<pos_t<Dim>>& shifts) {

    std::vector<pos_t<Dim>> res;
    res.reserve(shifts.size());
    for (auto shift : shifts) {
        auto new_pos = pos + shift;
        res.push_back(new_pos);
    }

    return res;
}


template <std::size_t Dim>
std::vector<pos_t<Dim>> make_poses(
    shifts_fn<Dim> shfn, const pos_t<Dim>& center, std::size_t range,
    std::optional<inside_fn<Dim>> infn) {

    return apply_shifts<Dim>(center, shfn(range), infn);
}


template <std::size_t Dim>
std::vector<pos_t<Dim>> make_poses(
    norm_fn<Dim> normfn, const pos_t<Dim>& center, std::size_t range,
    std::optional<inside_fn<Dim>> infn) {

    return apply_shifts<Dim>(center, make_shifts<Dim>(normfn, range), infn);
}


template <std::size_t Dim>
offsets_t make_offsets(
    shifts_fn<Dim> shfn, std::size_t center, const spt::vecu<Dim>& dimlens,
    std::size_t range, std::optional<inside_fn<Dim>> infn) {

    pos_t<Dim> center_pos = cgr::upos(center, dimlens);
    auto nbhpos = make_poses<Dim>(shfn, center_pos, range, infn);
    return poses_to_offsets<Dim>(nbhpos, dimlens);
}


template <std::size_t Dim>
offsets_t make_offsets(
    norm_fn<Dim> normfn, std::size_t center, const spt::vecu<Dim>& dimlens,
    std::size_t range, std::optional<inside_fn<Dim>> infn) {

    pos_t<Dim> center_pos = cgr::upos(center, dimlens);
    auto nbhpos = make_poses<Dim>(normfn, center_pos, range, infn);
    return poses_to_offsets<Dim>(nbhpos, dimlens);
}


template <std::size_t Dim, typename Cell>
using nbhood_t = std::vector<Cell*>;


template <std::size_t Dim, typename Cell>
nbhood_t<Dim, Cell> make_nbhood(
    norm_fn<Dim> normfn, const pos_t<Dim>& center, std::size_t range,
    try_cell_fn<Dim, Cell> trygetcell) {

    nbhood_t<Dim, Cell> res;
    auto nbhpos = trygetcell ? 
        make_poses<Dim>(normfn, center, range,
            [trygetcell](const pos_t<Dim>& pos) -> bool {
                return trygetcell.value()(pos); }) 
        : make_poses<Dim>(normfn, center, range);

    for (auto& pos : nbhpos) {
        auto pcell = trygetcell(pos);
        if (pcell)
            res.push_back(pcell);
    }

    return res;
}


template <std::size_t Dim, typename Cell>
nbhood_t<Dim, Cell> make_nbhood(
    const offsets_t& nbhoff,
    cell_by_offset_fn<Dim, Cell> getcell) {

    nbhood_t<Dim, Cell> res;
    res.reserve(nbhoff.size());
    for (auto& off : nbhoff)
        res.push_back(getcell(off));
    return res;
}


template <std::size_t Dim, typename Cell>
nbhood_t<Dim, Cell> make_nbhood(
    const std::vector<pos_t<Dim>>& nbhpos,
    try_cell_fn<Dim, Cell> trygetcell) {
    nbhood_t<Dim, Cell> res;
    for (auto& pos : nbhpos) {
        auto pcell = trygetcell(pos);
        if (pcell)
            res.push_back(pcell);
    }

    return res;
}

} // namespace cgr
