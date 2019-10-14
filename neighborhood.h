// Copyright © 2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <vector>
#include <functional>
#include "vec.h"


namespace cgr {

enum class neighborhood_type {
    von_neumann,
    moore
};


template <neighborhood_type NhoodType, std::size_t Dim>
class neighborhood_pos_impl;


template <std::size_t Dim>
class neighborhood_pos_impl<neighborhood_type::moore, Dim> {
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


template <std::size_t Dim>
class neighborhood_pos_impl<neighborhood_type::von_neumann, Dim> {
public:
    using veci = spt::veci<Dim>;
    using get_cell_func = std::function<Cell * (const veci&)>;

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


template <neighborhood_type NhoodType, std::size_t Dim, typename Cell>
class neighborhood {
public:
    using veci = spt::veci<Dim>;
    using get_cell_func = std::function<Cell*(const veci&)>;

    auto begin() const {
        return m_neighbors.cbegin();
    }
    auto end() const {
        return m_neighbors.cend();
    }
    Cell* central_cell() const {
        return m_central_cell;
    }

    neighborhood(const veci& center, std::size_t range, get_cell_func getcell) {
        m_central_cell = getcell(center);
        using nbhood_pos = neighborhood_pos_impl<NhoodType, Dim>;
        for (auto& pos : nbhood_pos::neighbors_pos(center)) {
            Cell* cell = getcell(pos);
            if (cell)
                m_neighbors.push_back(cell);
        }
    }


private:
    Cell* m_central_cell;
    std::vector<Cell*> m_neighbors;
};

} // namespace cgr
