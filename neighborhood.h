// Copyright � 2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <vector>
#include <functional>
#include "vec.h"


namespace cgr {

enum class nbhood_kind {
    von_neumann,
    moore
};


template <nbhood_kind NbhoodKind, std::size_t Dim, typename Cell>
class nbhood_pos_impl;


template <std::size_t Dim, typename Cell>
class nbhood_pos_impl<nbhood_kind::moore, Dim, Cell> {
public:
    using veci = spt::veci<Dim>;
    using get_cell_func = std::function<Cell*(const veci&)>;

    static std::vector<veci> neighbors_pos(const veci& center, std::size_t range) {
        std::vector<veci> res;
        if constexpr (Dim == 2) {
            res.reserve(4 * range * (range + 1)); // (2 * range + 1)^2 - 1
            for (std::int64_t y = -range; y <= range; y++) 
                for (std::int64_t x = -range; x <= range; x++) 
                    if (!(x == 0 && y == 0))
                        res.push_back(center + veci{x, y});

        } else if constexpr (Dim == 3) {
            static_assert(false);

        } else {
            static_assert(false);
        }
        
        return res;
    }
};


template <std::size_t Dim, typename Cell>
class nbhood_pos_impl<nbhood_kind::von_neumann, Dim, Cell> {
public:
    using veci = spt::veci<Dim>;
    using get_cell_func = std::function<Cell*(const veci&)>;

    static std::vector<veci> neighbors_pos(const veci& center, std::size_t range) {
        std::vector<veci> res;
        if constexpr (Dim == 2) {
            res.reserve(4 * range);
            for (std::int64_t y = -range; y <= range; y++)
                for (std::int64_t x = -range; x <= range; x++)
                    if (!(x == 0 && y == 0) &&
                        std::abs(x) + std::abs(y) <= range)
                        res.push_back(center + veci{x, y});

        } else if constexpr (Dim == 3) {
            static_assert(false);

        } else {
            static_assert(false);
        }

        return res;
    }
};


template <std::size_t Dim, typename Cell>
class nbhood {
public:
    using veci = spt::veci<Dim>;
    using get_cell_func = std::function<Cell*(const veci&)>;
    using get_pos_func = std::function<veci(const Cell*)>;

    auto begin() const {
        return m_neighbors.cbegin();
    }
    auto end() const {
        return m_neighbors.cend();
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
        nbhood_pos_impl<NbhoodKind, Dim, Cell>;
        std::vector<veci> neighbors_pos;
        switch (kind) {
        case nbhood_kind::von_neumann:
            neighbors_pos = nbhood_pos_impl<nbhood_kind::von_neumann, 
                Dim, Cell>::neighbors_pos(center);
            break;
        case nbhood_kind::moore:
            neighbors_pos = nbhood_pos_impl<nbhood_kind::moore,
                Dim, Cell>::neighbors_pos(center);
            break;
        }
        for (auto& pos : neighbors_pos) {
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
