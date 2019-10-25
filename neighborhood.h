// Copyright © 2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <vector>
#include <functional>
#include "vec.h"


namespace cgr {

enum class nbhood_kind {
    von_neumann,
    moore,
    euclid
};


template <nbhood_kind NbhoodKind, std::size_t Dim, 
    typename TryGetCell = std::function<bool(const spt::veci<Dim>&)>>
class nbhood_pos_impl;


template <std::size_t Dim, typename TryGetCell>
class nbhood_pos_impl<nbhood_kind::moore, Dim, TryGetCell> {
public:
    using veci = spt::veci<Dim>;
    static std::vector<veci> run(const veci& center, std::size_t range,
                                 std::optional<TryGetCell> getcell = std::nullopt) {
        std::int64_t srange = range;
        std::vector<veci> res;
        if constexpr (Dim == 2) {
            res.reserve(4 * range * (range + 1)); // (2 * range + 1)^2 - 1
            for (std::int64_t y = -srange; y <= srange; ++y)
                for (std::int64_t x = -srange; x <= srange; ++x)
                    if (!(x == 0 && y == 0)) {
                        auto new_pos = center + veci{x, y};
                        if (!getcell) 
                            res.push_back(new_pos);
                        else
                            if (getcell.value()(new_pos))
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
                            if (!getcell)
                                res.push_back(new_pos);
                            else
                                if (getcell.value()(new_pos))
                                    res.push_back(new_pos);
                        }

        } else {
            static_assert(false);
        }
        
        return res;
    }
};


template <std::size_t Dim, typename TryGetCell>
class nbhood_pos_impl<nbhood_kind::von_neumann, Dim, TryGetCell> {
public:
    using veci = spt::veci<Dim>;
    static std::vector<veci> run(const veci& center, std::size_t range,
                                 std::optional<TryGetCell> getcell = std::nullopt) {
        std::int64_t srange = range;
        std::vector<veci> res;
        if constexpr (Dim == 2) {
            res.reserve(4 * range);
            for (std::int64_t y = -srange; y <= srange; ++y)
                for (std::int64_t x = -srange; x <= srange; ++x)
                    if (!(x == 0 && y == 0) &&
                        std::abs(x) + std::abs(y) <= srange) {
                        auto new_pos = center + veci{x, y};
                        if (!getcell) 
                            res.push_back(new_pos);
                        else 
                            if (getcell.value()(new_pos))
                                res.push_back(new_pos);
                    }

        } else if constexpr (Dim == 3) {
            std::size_t buf = 2 * range + 1;
            res.reserve(buf * buf * buf - 1);
            for (std::int64_t z = -srange; z <= srange; ++z)
                for (std::int64_t y = -srange; y <= srange; ++y)
                    for (std::int64_t x = -srange; x <= srange; ++x)
                        if (!(x == 0 && y == 0 && z == 0) &&
                            std::abs(x) + std::abs(y) + std::abs(z) <= srange) {
                            auto new_pos = center + veci{x, y, z};
                            if (!getcell)
                                res.push_back(new_pos);
                            else
                                if (getcell.value()(new_pos))
                                    res.push_back(new_pos);
                        }

        } else {
            static_assert(false);
        }

        return res;
    }
};


template <std::size_t Dim, typename TryGetCell>
class nbhood_pos_impl<nbhood_kind::euclid, Dim, TryGetCell> {
public:
    using veci = spt::veci<Dim>;
    static std::vector<veci> run(const veci& center, std::size_t range, 
                                 std::optional<TryGetCell> getcell = std::nullopt) {
        std::int64_t srange = range;
        std::int64_t srange2 = srange * srange;
        std::vector<veci> res;
        if constexpr (Dim == 2) {
            for (std::int64_t y = -srange; y <= srange; ++y)
                for (std::int64_t x = -srange; x <= srange; ++x)
                    if (!(x == 0 && y == 0) &&
                        x * x + y * y <= srange2) {
                        auto new_pos = center + veci{x, y};
                        if (!getcell)
                            res.push_back(new_pos);
                        else
                            if (getcell.value()(new_pos))
                                res.push_back(new_pos);
                    }

        } else if constexpr (Dim == 3) {
            for (std::int64_t z = -srange; z <= srange; ++z)
                for (std::int64_t y = -srange; y <= srange; ++y)
                    for (std::int64_t x = -srange; x <= srange; ++x)
                        if (!(x == 0 && y == 0 && z == 0) &&
                            x * x + y * y + z * z <= srange2) {
                            auto new_pos = center + veci{x, y, z};
                            if (!getcell)
                                res.push_back(new_pos);
                            else
                                if (getcell.value()(new_pos))
                                    res.push_back(new_pos);
                        }

        } else {
            static_assert(false);
        }

        return res;
    }
};


template <std::size_t Dim>
using nbhood_pos = std::vector<spt::veci<Dim>>;


template <std::size_t Dim>
nbhood_pos<Dim> make_nbhood_pos(const spt::veci<Dim>& center, nbhood_kind kind, std::size_t range, 
                                std::optional<std::function<bool(const spt::veci<Dim>&)>> getcell = std::nullopt) {
    switch (kind) {
    case nbhood_kind::von_neumann:
        return nbhood_pos_impl<nbhood_kind::von_neumann, Dim>
            ::run(center, range, getcell);

    case nbhood_kind::moore:
        return nbhood_pos_impl<nbhood_kind::moore, Dim>
            ::run(center, range, getcell);

    case nbhood_kind::euclid:
        return nbhood_pos_impl<nbhood_kind::euclid, Dim>
            ::run(center, range, getcell);

    default:
        return {};
    }
}


template <std::size_t Dim, typename Cell>
using nbhood = std::vector<Cell*>;


template <std::size_t Dim, typename Cell>
nbhood<Dim, Cell> make_nbhood(const spt::veci<Dim>& center, nbhood_kind kind, std::size_t range, 
                              std::function<Cell*(const spt::veci<Dim>&)> getcell) {
    nbhood<Dim, Cell> res;
    auto nbhpos = getcell ? make_nbhood_pos(center, kind, range, 
                                            [getcell](const spt::veci<Dim>& pos) 
                                            -> bool { return getcell.value()(pos); }) 
        : make_nbhood_pos(center, kind, range);

    for (auto& nbpos : nbhpos) {
        auto pcell = getcell(nbpos);
        if (pcell)
            res.push_back(pcell);
    }

    return res;
}


template <std::size_t Dim, typename Cell>
    nbhood<Dim, Cell> make_nbhood(const nbhood_pos<Dim>& nbhpos,
                                  std::function<Cell*(const spt::veci<Dim>&)> getcell) {
    nbhood<Dim, Cell> res;
    for (auto& nbpos : nbhpos) {
        auto pcell = getcell(nbpos);
        if (pcell)
            res.push_back(pcell);
    }

    return res;
}

} // namespace cgr
