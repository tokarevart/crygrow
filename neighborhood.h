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


template <nbhood_kind NbhoodKind, std::size_t Dim>
class nbhood_pos_impl;


template <std::size_t Dim>
class nbhood_pos_impl<nbhood_kind::moore, Dim> {
public:
    using veci = spt::veci<Dim>;

    static std::vector<veci> neighbors_pos(const veci& center, std::size_t range) {
        std::int64_t srange = range;
        std::vector<veci> res;
        if constexpr (Dim == 2) {
            res.reserve(4 * range * (range + 1)); // (2 * range + 1)^2 - 1
            for (std::int64_t y = -srange; y <= srange; y++)
                for (std::int64_t x = -srange; x <= srange; x++)
                    if (!(x == 0 && y == 0))
                        res.push_back(center + veci{x, y});

        } else if constexpr (Dim == 3) {
            std::size_t buf = 2 * range + 1;
            res.reserve(buf * buf * buf - 1);
            for (std::int64_t z = -srange; z <= srange; z++)
                for (std::int64_t y = -srange; y <= srange; y++)
                    for (std::int64_t x = -srange; x <= srange; x++)
                        if (!(x == 0 && y == 0 && z == 0))
                            res.push_back(center + veci{x, y, z});

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

    static std::vector<veci> neighbors_pos(const veci& center, std::size_t range) {
        std::int64_t srange = range;
        std::vector<veci> res;
        if constexpr (Dim == 2) {
            res.reserve(4 * range);
            for (std::int64_t y = -srange; y <= srange; y++)
                for (std::int64_t x = -srange; x <= srange; x++)
                    if (!(x == 0 && y == 0) &&
                        std::abs(x) + std::abs(y) <= srange)
                        res.push_back(center + veci{x, y});

        } else if constexpr (Dim == 3) {
            std::size_t buf = 2 * range + 1;
            res.reserve(buf * buf * buf - 1);
            for (std::int64_t z = -srange; z <= srange; z++)
                for (std::int64_t y = -srange; y <= srange; y++)
                    for (std::int64_t x = -srange; x <= srange; x++)
                        if (!(x == 0 && y == 0 && z == 0) &&
                            std::abs(x) + std::abs(y) + std::abs(z) <= srange)
                            res.push_back(center + veci{x, y, z});

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

    static std::vector<veci> neighbors_pos(const veci& center, std::size_t range) {
        std::int64_t srange = range;
        std::vector<veci> res;
        if constexpr (Dim == 2) {
            for (std::int64_t y = -srange; y <= srange; y++)
                for (std::int64_t x = -srange; x <= srange; x++)
                    if (!(x == 0 && y == 0) &&
                        x * x + y * y <= srange * srange)
                        res.push_back(center + veci{ x, y });

        } else if constexpr (Dim == 3) {
            for (std::int64_t z = -srange; z <= srange; z++)
                for (std::int64_t y = -srange; y <= srange; y++)
                    for (std::int64_t x = -srange; x <= srange; x++)
                        if (!(x == 0 && y == 0 && z == 0) &&
                            x * x + y * y + z * z <= srange * srange)
                            res.push_back(center + veci{ x, y, z });

        } else {
            static_assert(false);
        }

        return res;
    }
};


template <std::size_t Dim>
std::vector<spt::veci<Dim>> neighbors_pos(const spt::veci<Dim>& center, nbhood_kind kind, std::size_t range) {
    switch (kind) {
    case nbhood_kind::von_neumann:
        return nbhood_pos_impl<nbhood_kind::von_neumann,
            Dim>::neighbors_pos(center, range);

    case nbhood_kind::moore:
        return nbhood_pos_impl<nbhood_kind::moore,
            Dim>::neighbors_pos(center, range);

    case nbhood_kind::euclid:
        return nbhood_pos_impl<nbhood_kind::euclid,
            Dim>::neighbors_pos(center, range);
    }
}


template <std::size_t Dim, typename Cell>
class nbhood {
public:
    using veci = spt::veci<Dim>;
    using get_pos_func = std::function<veci(const Cell*)>;
    using get_cell_func = std::function<Cell*(const veci&)>;

    auto begin() const {
        return m_neighbors.cbegin();
    }
    auto end() const {
        return m_neighbors.cend();
    }
    bool empty() const {
        return m_neighbors.empty();
    }
    std::size_t size() const {
        return m_neighbors.size();
    }
    Cell* central_cell() const {
        return m_central_cell;
    }
    nbhood_kind kind() const {
        return m_kind;
    }
    std::size_t range() const {
        return m_range;
    }

    nbhood(const veci& center, nbhood_kind kind, std::size_t range, get_cell_func getcell)
        : m_central_cell{getcell(center)}, m_kind{kind}, m_range{range} {
        for (auto& pos : neighbors_pos(center, kind, range)) {
            Cell* cell = getcell(pos);
            if (cell)
                m_neighbors.push_back(cell);
        }
    }
    nbhood(const Cell* central_cell, nbhood_kind kind, std::size_t range,
                 get_pos_func getpos, get_cell_func getcell)
        : nbhood(getpos(central_cell), kind, range, getcell) {}


private:
    Cell* m_central_cell;
    nbhood_kind m_kind;
    std::size_t m_range;
    std::vector<Cell*> m_neighbors;
};

} // namespace cgr
