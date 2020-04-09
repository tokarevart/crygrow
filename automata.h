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
#include "clr-grain.h"


namespace cgr {

template <std::size_t Dim, typename Real = double>
class automata {
public:
    static constexpr std::size_t dim = Dim;
    using cell_type = cgr::cell<Dim, Real>;
    using grain_type = typename cell_type::grain_type;
    using material_type = typename grain_type::material_type;
    using orientation_type = typename grain_type::orientation_type;
    using clr_grain_type = clr_grain<Dim, Real>;
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
    const std::vector<cell_type>& cells() const {
        return m_cells;
    }

    const cell_type& cell(std::size_t offset) const {
        return m_cells[offset];
    }
    const cell_type& cell(const upos_t<Dim>& pos) const {
        return cell(offset(pos));
    }
    const cell_type& cell(const pos_t<Dim>& pos) const {
        return cell(offset(pos));
    }
    void set_cell(std::size_t offset, const cell_type& new_cell) {
        m_cells[offset] = new_cell;
    }
    void set_cell(const pos_t<Dim>& pos, const cell_type& new_cell) {
        m_cells[offset(pos)] = new_cell;
    }
    template <typename PosesContainer, typename CellsContainer>
    void set_cells(const PosesContainer& poses, const CellsContainer& cells) {
        auto pos_it = poses.begin();
        auto cell_it = cells.begin();
        for (; pos_it != poses.end(); ++pos_it, ++cell_it)
            set_cell(*pos_it, *cell_it);
    }

    bool try_set_cell(const pos_t<Dim>& pos, const cell_type& new_cell) {
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
        for (auto& clrg : m_clrgrains)
            clrg.set_range(m_range);
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

        for (auto& clrg : m_clrgrains)
            for (std::size_t off : clrg.front())
                m_cells[off].crysted = true;

        #pragma omp parallel for
        for (std::int64_t i = 0; i < m_clrgrains.size(); ++i)
            m_clrgrains[i].advance_front(
                [this](std::size_t off) -> bool { return m_cells[off].crysted; },
                [this](std::size_t off) -> std::size_t { return m_cells[off].grains.size(); });

        for (auto& clrg : m_clrgrains)
            for (std::size_t off : clrg.front())
                m_cells[off].grains.push_back(clrg.grain());

        return true;
    }

    void spawn_grain(const grain_type* grain, const upos_t<Dim>& nucleus_pos, nbh::nbhood_kind kind) {
        spawn_grain(grain, offset(nucleus_pos), kind);
    }
    void spawn_grain(const grain_type* grain, std::size_t nucleus_off, nbh::nbhood_kind kind) {
        m_clrgrains.emplace_back(grain, kind, m_dim_lens, nucleus_off);
        m_clrgrains.back().set_range(m_range);
        m_cells[nucleus_off].crysted = true;
        m_cells[nucleus_off].grains.push_back(grain);
    }

    automata(std::size_t dimlen)
        : automata(upos_t<Dim>::filled_with(dimlen)) {}
    automata(const upos_t<Dim>& dimlens) : m_dim_lens{ dimlens } {
        std::size_t new_num_cells = std::accumulate(
            m_dim_lens.x.begin(), m_dim_lens.x.end(),
            static_cast<std::size_t>(1), std::multiplies<std::size_t>());
        m_cells.assign(new_num_cells, cell_type());
    }


private:
    std::size_t m_range = 0;
    upos_t<Dim> m_dim_lens;

    std::vector<cell_type> m_cells;
    std::vector<clr_grain_type> m_clrgrains;
    // todo: store common cells, such as empty cell and crystallized cell for each grain
};

} // namespace cgr
