// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <algorithm>
#include <map>
#include <cstddef>
#include <memory>
#include <optional>
#include "neighborhood.h"
#include "vec.h"
#include "sptops.h"
#include "cgralgs.h"
#include "cell.h"


namespace cgr {

template <std::size_t Dim, typename Real = double>
class automata {
public:
    static constexpr std::size_t dim = Dim;
    using cell_type = cgr::cell<Dim, Real>;
    using cells_container = std::vector<cell_type>;
    using grain_type = typename cell_type::grain_type;
    using material_type = typename grain_type::material_type;
    using orientation_type = typename grain_type::orientation_type;
    using grow_dir_type = grow_dir_t<Dim, Real>;

    std::size_t num_crysted_cells() const {
        std::size_t res = 0;
        for (std::size_t i = 0; i < num_cells(); ++i)
            if (cell(i).crysted)
                ++res;
        return res;
    }
    std::size_t num_cells() const {
        return m_cells.size();
    }
    const cells_container& cells() const {
        return m_cells;
    }

    cell_type cell(std::size_t offset) const {
        return m_cells[offset];
    }
    cell_type cell(const upos_t<Dim>& pos) const {
        return cell(offset(pos));
    }
    cell_type cell(const pos_t<Dim>& pos) const {
        return cell(offset(pos));
    }
    void set_cell(std::size_t offset, cell_type new_cell) {
        m_cells[offset] = new_cell;
    }
    void set_cell(const pos_t<Dim>& pos, cell_type new_cell) {
        m_cells[offset(pos)] = new_cell;
    }
    template <typename PosesContainer, typename CellsContainer>
    void set_cells(const PosesContainer& poses, const CellsContainer& cells) {
        auto pos_it = poses.begin();
        auto cell_it = cells.begin();
        for (; pos_it != poses.end(); ++pos_it, ++cell_it)
            set_cell(*pos_it, *cell_it);
    }

    std::optional<cell_type> try_cell(const upos_t<Dim>& pos) const {
        return inside(pos) ? cell(pos) : std::nullopt;
    }
    std::optional<cell_type> try_cell(const pos_t<Dim>& pos) const {
        return inside(pos) ? cell(pos) : std::nullopt;
    }
    bool try_set_cell(const pos_t<Dim>& pos, cell_type new_cell) {
        if (inside(pos)) {
            m_cells[offset(pos)] = new_cell;
            return true;
        }

        return false;
    }
    template <typename PosesContainer, typename CellsContainer>
    void try_set_cells(const PosesContainer& poses, const CellsContainer& cells) {
        auto pos_it = poses.begin();
        auto cell_it = cells.begin();
        for (; pos_it != poses.end(); ++pos_it, ++cell_it)
            try_set_cell(*pos_it, *cell_it);
    }

    bool inside(const upos_t<Dim>& pos) const {
        return cgr::inside(pos, m_dim_lens);
    }
    bool inside(const pos_t<Dim>& pos) const {
        for (auto& e : pos.x)
            if (e < 0)
                return false;
        return cgr::inside(static_cast<upos_t<Dim>>(pos), m_dim_lens);
    }

    upos_t<Dim> upos(std::size_t offset) const {
        return cgr::upos(offset, m_dim_lens);
    }
    std::size_t offset(const upos_t<Dim>& pos) const {
        return offset(static_cast<pos_t<Dim>>(pos));
    }
    std::size_t offset(const pos_t<Dim>& pos) const {
        return cgr::offset(pos, m_dim_lens);
    }

    const upos_t<Dim>& dim_lens() const {
        return m_dim_lens;
    }
    std::size_t range() const {
        return m_range;
    }
    void set_range(std::size_t range) {
        m_range = range;
    }

    bool stop_condition() const {
        for (std::size_t i = 0; i < num_cells(); ++i)
            if (!cell(i).crysted)
                return false;
        return true;
    }
    bool iterate() {
        if (stop_condition())
            return false;

        //

        return true;
    }

    automata(std::size_t dimlen)
        : automata(upos_t<Dim>::filled_with(dimlen)) {}
    automata(const upos_t<Dim>& dimlens) {
        m_dim_lens = dimlens;
        std::size_t new_num_cells = std::accumulate(
            m_dim_lens.x.begin(), m_dim_lens.x.end(),
            static_cast<std::size_t>(1), std::multiplies<std::size_t>());
        m_cells.assign(new_num_cells, cell_type());
    }


private:
    std::size_t m_range = 0;
    upos_t<Dim> m_dim_lens;

    cells_container m_cells;
    // todo: store common cells, such as empty cell and crystallized cell for each grain
};

} // namespace cgr
