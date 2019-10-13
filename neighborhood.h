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


template <neighborhood_type NhoodType, std::size_t Dim, std::size_t Range, typename Cell>
class neighborhood_impl {

};


template <neighborhood_type NhoodType, std::size_t Dim, std::size_t Range, typename Cell>
class neighborhood {
public:
    using veci = spt::veci<Dim>;
    using cell_func = std::function<Cell*(const veci&)>;

    auto begin() const {
        return m_neighbors.cbegin();
    }
    auto end() const {
        return m_neighbors.cend();
    }

    neighborhood(const veci& center, cell_func cellf) {

    }


private:
    std::vector<Cell*> m_neighbors;
};

} // namespace cgr
