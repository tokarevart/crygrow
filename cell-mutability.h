// Copyright © 2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once

namespace cgr {

enum class cell_mut_group {
    mutable_only,
    universal
};

enum class cell_mut {
    mutable_cell,
    constant_cell
};


template <cell_mut_group CellMut>
class cell_base;


template <>
class cell_base<cell_mut_group::mutable_only> {
public:
    static constexpr cell_mut_group mutability_group = cell_mut_group::mutable_only;
};


template <>
class cell_base<cell_mut_group::universal> {
public:
    static constexpr cell_mut_group mutability_group = cell_mut_group::universal;

    cell_mut mutability() const {
        return m_mut;
    }

    cell_base(cell_mut mut) : m_mut{mut} {}


private:
    cell_mut m_mut;
};

} // namespace cgr
