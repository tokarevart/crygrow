// Copyright � 2019 Tokarev Artem. All rights reserved.
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

template <nbhood_kind NbhoodKind, std::size_t Dim>
bool inside_nbhood(const spt::veci<Dim>& pos, std::size_t range) {
    if constexpr (NbhoodKind == nbhood_kind::von_neumann) {
        return cgr::magnitude_von_neumann(pos) <= range;

    } else if constexpr (NbhoodKind == nbhood_kind::moore) {
        for (auto e : pos.x)
            if (e > static_cast<std::int64_t>(range))
                return false;
        return true;

    } else if constexpr (NbhoodKind == nbhood_kind::euclid) {
        return cgr::magnitude2_euclid(pos) <= range * range;
    }
}


using nbhood_offset = std::vector<std::size_t>;

template <std::size_t Dim, typename Cell>
using try_get_cell_t = std::function<Cell*(const spt::veci<Dim>&)>;
template <std::size_t Dim>
using try_get_bcell_t = std::function<bool(const spt::veci<Dim>&)>;


template <std::size_t Dim>
using nbhood_pos = std::vector<spt::veci<Dim>>;


template <nbhood_kind NbhoodKind, std::size_t Dim>
class nbhood_pos_impl;


template <std::size_t Dim>
class nbhood_pos_impl<nbhood_kind::moore, Dim> {
public:
    using veci = spt::veci<Dim>;
    static nbhood_pos<Dim> run(const veci& center, std::size_t range,
                               std::optional<try_get_bcell_t<Dim>> trygetcell) {
        std::int64_t srange = range;
        std::vector<veci> res;
        if constexpr (Dim == 2) {
            res.reserve(4 * range * (range + 1)); // (2 * range + 1)^2 - 1
            for (std::int64_t y = -srange; y <= srange; ++y)
                for (std::int64_t x = -srange; x <= srange; ++x)
                    if (!(x == 0 && y == 0)) {
                        auto new_pos = center + veci{x, y};
                        if (!trygetcell) 
                            res.push_back(new_pos);
                        else
                            if (trygetcell.value()(new_pos))
                                res.push_back(new_pos);
                    }

        } else if constexpr (Dim == 3) {
            std::size_t buf = 2 * range + 1;
            res.reserve(buf * buf * buf - 1);
            for (std::int64_t z = -srange; z <= srange; ++z)
                for (std::int64_t y = -srange; y <= srange; ++y)
                    for (std::int64_t x = -srange; x <= srange; ++x)
                        if (!(x == 0 && y == 0 && z == 0)) {
                            auto new_pos = center + veci{x, y, z};
                            if (!trygetcell)
                                res.push_back(new_pos);
                            else
                                if (trygetcell.value()(new_pos))
                                    res.push_back(new_pos);
                        }

        } else {
            static_assert(false);
        }
        
        return res;
    }
};


template <std::size_t Dim>
class nbhood_pos_impl<nbhood_kind::von_neumann, Dim> {
public:
    using veci = spt::veci<Dim>;
    static nbhood_pos<Dim> run(const veci& center, std::size_t range,
                               std::optional<try_get_bcell_t<Dim>> trygetcell) {
        std::int64_t srange = range;
        std::vector<veci> res;
        if constexpr (Dim == 2) {
            res.reserve(4 * range * (range + 1)); // (2 * range + 1)^2 - 1
            for (std::int64_t y = -srange; y <= srange; ++y)
                for (std::int64_t x = -srange; x <= srange; ++x)
                    if (!(x == 0 && y == 0) &&
                        std::abs(x) + std::abs(y) <= srange) {
                        auto new_pos = center + veci{x, y};
                        if (!trygetcell) 
                            res.push_back(new_pos);
                        else 
                            if (trygetcell.value()(new_pos))
                                res.push_back(new_pos);
                    }
            res.shrink_to_fit();

        } else if constexpr (Dim == 3) {
            std::size_t buf = 2 * range + 1;
            res.reserve(buf * buf * buf - 1);
            for (std::int64_t z = -srange; z <= srange; ++z)
                for (std::int64_t y = -srange; y <= srange; ++y)
                    for (std::int64_t x = -srange; x <= srange; ++x)
                        if (!(x == 0 && y == 0 && z == 0) &&
                            std::abs(x) + std::abs(y) + std::abs(z) <= srange) {
                            auto new_pos = center + veci{x, y, z};
                            if (!trygetcell)
                                res.push_back(new_pos);
                            else
                                if (trygetcell.value()(new_pos))
                                    res.push_back(new_pos);
                        }
            res.shrink_to_fit();

        } else {
            static_assert(false);
        }

        return res;
    }
};


template <std::size_t Dim>
class nbhood_pos_impl<nbhood_kind::euclid, Dim> {
public:
    using veci = spt::veci<Dim>;
    static nbhood_pos<Dim> run(const veci& center, std::size_t range,
                               std::optional<try_get_bcell_t<Dim>> trygetcell) {
        std::int64_t srange = range;
        std::int64_t srange2 = srange * srange;
        std::vector<veci> res;
        if constexpr (Dim == 2) {
            res.reserve(4 * range * (range + 1)); // (2 * range + 1)^2 - 1
            for (std::int64_t y = -srange; y <= srange; ++y)
                for (std::int64_t x = -srange; x <= srange; ++x)
                    if (!(x == 0 && y == 0) &&
                        x * x + y * y <= srange2) {
                        auto new_pos = center + veci{x, y};
                        if (!trygetcell)
                            res.push_back(new_pos);
                        else
                            if (trygetcell.value()(new_pos))
                                res.push_back(new_pos);
                    }
            res.shrink_to_fit();

        } else if constexpr (Dim == 3) {
            std::size_t buf = 2 * range + 1;
            res.reserve(buf * buf * buf - 1);
            for (std::int64_t z = -srange; z <= srange; ++z)
                for (std::int64_t y = -srange; y <= srange; ++y)
                    for (std::int64_t x = -srange; x <= srange; ++x)
                        if (!(x == 0 && y == 0 && z == 0) &&
                            x * x + y * y + z * z <= srange2) {
                            auto new_pos = center + veci{x, y, z};
                            if (!trygetcell)
                                res.push_back(new_pos);
                            else
                                if (trygetcell.value()(new_pos))
                                    res.push_back(new_pos);
                        }
            res.shrink_to_fit();

        } else {
            static_assert(false);
        }

        return res;
    }
};


template <nbhood_kind NbhoodKind, std::size_t Dim>
nbhood_offset make_nbhood_offset(std::size_t center, const spt::vecu<Dim>& dim_lens, std::size_t range,
                                 std::optional<try_get_bcell_t<Dim>> trygetcell = std::nullopt) {
    spt::veci<Dim> center_upos = cgr::upos(center, dim_lens);
    nbhood_offset res;

    auto nbdimlens = spt::vecu<Dim>::filled_with(2 * range + 1);
    std::size_t num_nboffsets = std::accumulate(nbdimlens.x.begin(), nbdimlens.x.end(),
                                                static_cast<std::size_t>(1),
                                                std::multiplies<std::size_t>());
    auto nbcenter = nbdimlens / 2;
    res.reserve(num_nboffsets - 1);
    for (std::size_t i = 0; i < num_nboffsets; ++i) {
        if (i == num_nboffsets / 2)
            continue;

        auto relpos = static_cast<spt::veci<Dim>>(cgr::upos(i, nbdimlens)) - nbcenter;
        if (inside_nbhood<NbhoodKind, Dim>(relpos, range)) {
            auto new_pos = center_upos + relpos;
            if (!trygetcell)
                res.push_back(cgr::offset(new_pos, dim_lens));
            else
                if (trygetcell.value()(new_pos))
                    res.push_back(cgr::offset(new_pos, dim_lens));
        }
    }
    res.shrink_to_fit();

    return res;
}


template <nbhood_kind NbhoodKind, std::size_t Dim>
nbhood_pos<Dim> make_nbhood_pos(const spt::veci<Dim>& center, std::size_t range, 
                                std::optional<try_get_bcell_t<Dim>> trygetcell = std::nullopt) {
    return nbhood_pos_impl<NbhoodKind, Dim>::run(center, range, trygetcell);
}


template <std::size_t Dim, typename Cell>
using nbhood = std::vector<Cell*>;


template <std::size_t Dim, typename Cell>
nbhood<Dim, Cell> make_nbhood(std::size_t center, const spt::vecu<Dim>& dim_lens, 
                              nbhood_kind kind, std::size_t range,
                              try_get_cell_t<Dim, Cell> trygetcell) {
    nbhood<Dim, Cell> res;
    auto nbhpos = trygetcell ? make_nbhood_pos(center, kind, range,
                                               [trygetcell](std::size_t pos)
                                               -> bool { return trygetcell.value()(pos); })
        : make_nbhood_pos(center, kind, range);

    for (auto& nbpos : nbhpos) {
        auto pcell = trygetcell(nbpos);
        if (pcell)
            res.push_back(pcell);
    }

    return res;
}


template <std::size_t Dim, typename Cell>
nbhood<Dim, Cell> make_nbhood(const spt::veci<Dim>& center, nbhood_kind kind, std::size_t range, 
                              try_get_cell_t<Dim, Cell> trygetcell) {
    nbhood<Dim, Cell> res;
    auto nbhpos = trygetcell ? make_nbhood_pos(center, kind, range, 
                                            [trygetcell](const spt::veci<Dim>& pos) 
                                            -> bool { return trygetcell.value()(pos); }) 
        : make_nbhood_pos(center, kind, range);

    for (auto& nbpos : nbhpos) {
        auto pcell = trygetcell(nbpos);
        if (pcell)
            res.push_back(pcell);
    }

    return res;
}


template <std::size_t Dim, typename Cell>
nbhood<Dim, Cell> make_nbhood(const nbhood_offset& nbhpos,
                              try_get_cell_t<Dim, Cell> trygetcell) {
    nbhood<Dim, Cell> res;
    for (auto& nbpos : nbhpos) {
        auto pcell = trygetcell(nbpos);
        if (pcell)
            res.push_back(pcell);
    }

    return res;
}


template <std::size_t Dim, typename Cell>
nbhood<Dim, Cell> make_nbhood(const nbhood_pos<Dim>& nbhpos,
                              try_get_cell_t<Dim, Cell> trygetcell) {
    nbhood<Dim, Cell> res;
    for (auto& nbpos : nbhpos) {
        auto pcell = trygetcell(nbpos);
        if (pcell)
            res.push_back(pcell);
    }

    return res;
}

} // namespace cgr
