// Copyright � 2019-2020 Artyom Tokarev. All rights reserved.
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
    static constexpr Real epsilon = std::numeric_limits<Real>::epsilon();
    using cell_type = cgr::cell<Dim, Real>;
    using cells_container = std::vector<cell_type*>;
    using nbhood_type = cgr::nbh::nbhood_t<Dim, cell_type>;
    using nbhood_pos_type = cgr::nbh::poses_t<Dim>;
    using nbhood_offset_type = cgr::nbh::offsets_t;
    using nbhood_offsets_container = std::vector<nbhood_offset_type>;
    using cells_delta_container = std::vector<Real>;
    using grain_type = typename cell_type::grain_type;
    using material_type = typename grain_type::material_type;
    using orientation_type = typename grain_type::orientation_type;
    using grow_dir_type = grow_dir_t<Dim, Real>;

    std::size_t num_crystallized_cells() const {
        std::size_t res = 0;
        for (std::size_t i = 0; i < num_cells(); ++i)
            if (cell(i)->crystallinity >= 1.0 - epsilon * (1.0 + cell(i)->crystallinity))
                ++res;
        return res;
    }
    std::size_t num_cells() const {
        return m_cells.size();
    }
    const cells_container& cells() const {
        return m_cells;
    }

    cell_type* cell(std::size_t offset) const {
        return m_cells[offset];
    }
    cell_type* cell(const upos_t<Dim>& pos) const {
        return cell(offset(pos));
    }
    cell_type* cell(const pos_t<Dim>& pos) const {
        return cell(offset(pos));
    }
    void set_cell(std::size_t offset, const cell_type* new_cell) {
        m_cells[offset] = const_cast<cell_type*>(new_cell);
    }
    void set_cell(const pos_t<Dim>& pos, const cell_type* new_cell) {
        m_cells[offset(pos)] = const_cast<cell_type*>(new_cell);
    }
    template <typename PosesContainer, typename CellsContainer>
    void set_cells(const PosesContainer& poses, const CellsContainer& cells) {
        auto pos_it = poses.begin();
        auto cell_it = cells.begin();
        for (; pos_it != poses.end(); ++pos_it, ++cell_it)
            set_cell(*pos_it, &(*cell_it));
    }

    cell_type* try_cell(const upos_t<Dim>& pos) const {
        return inside(pos) ? cell(pos) : nullptr;
    }
    cell_type* try_cell(const pos_t<Dim>& pos) const {
        return inside(pos) ? cell(pos) : nullptr;
    }
    bool try_set_cell(const pos_t<Dim>& pos, const cell_type* new_cell) {
        if (inside(pos)) {
            m_cells[offset(pos)] = const_cast<cell_type*>(new_cell);
            return true;
        }

        return false;
    }
    template <typename PosesContainer, typename CellsContainer>
    void try_set_cells(const PosesContainer& poses, const CellsContainer& cells) {
        auto pos_it = poses.begin();
        auto cell_it = cells.begin();
        for (; pos_it != poses.end(); ++pos_it, ++cell_it)
            try_set_cell(*pos_it, &(*cell_it));
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
        return cgr::offset(static_cast<pos_t<Dim>>(pos), m_dim_lens);
    }

    nbhood_offset_type make_offsets(std::size_t pos_offset) const {
        return cgr::make_offsets<Dim>(m_default_shift_fn,
            pos_offset, m_dim_lens, m_default_range,
            [this](const pos_t<Dim>& pos) -> bool { return inside(pos); });
    }
    const nbhood_offset_type& get_nbhood_offset(std::size_t pos_offset) const {
        return m_nbhood_offsets[pos_offset];
    }
    void set_nbhood_offset(std::size_t pos_offset) {
        m_nbhood_offsets[pos_offset] = make_offsets(pos_offset);
    }
    void set_nbhood_offset(std::size_t pos_offset, nbhood_offset_type&& nbh_offset) {
        m_nbhood_offsets[pos_offset] = std::move(nbh_offset);
    }
    void clear_nbhood_offset(std::size_t pos_offset) {
        m_nbhood_offsets[pos_offset].clear();
        m_nbhood_offsets[pos_offset].shrink_to_fit();
    }

    void reserve(std::size_t count) {
        reserve_cells(count);
        reserve_nbhood_offsets(count);
    }
    void reserve_cells(std::size_t count) {
        m_cells.reserve(count);
    }
    void reserve_nbhood_offsets(std::size_t count) {
        m_nbhood_offsets.reserve(count);
    }

    void erase(const pos_t<Dim>& pos) {
        erase_cell(pos);
        erase_nbhood_poses(pos);
    }
    void erase_cell(const pos_t<Dim>& pos) {
        set_cell(pos, nullptr);
    }
    void erase_nbhood_poses(const pos_t<Dim>& pos) {
        if (inside(pos))
            m_nbhood_poses[offset(pos)].clear();
    }

    const upos_t<Dim>& dim_lens() const {
        return m_dim_lens;
    }

    std::size_t default_range() const {
        return m_default_range;
    }
    std::size_t default_nbhood_size() const {
        return m_default_nbhood_size;
    }

    void set_default_range(std::size_t range) {
        m_default_range = range;
        set_default_nbhood_size();
    }

    bool stop_condition() const {
        for (std::size_t i = 0; i < num_cells(); ++i)
            if (cell(i)->crystallinity < 1.0 - epsilon * (1.0 + cell(i)->crystallinity))
                return false;

        return true;
    }
    bool iterate() {
        if (stop_condition())
            return false;

        if (!is_nbhood_offsets_initialized)
            initialize_nbhood_offsets();

        for (auto& delta : m_cells_delta)
            delta = 0.0;

        auto defnbhoffset = make_offsets();
        std::int64_t accdefdpmagn2 = 0;
        for (auto nboff : defnbhoffset)
            accdefdpmagn2 += (
                pos_t<Dim>::filled_with(default_range())
                - cgr::upos(nboff, upos_t<Dim>::filled_with(2 * default_range() + 1))
                ).magnitude2();
        
        #pragma omp parallel 
        {
            #pragma omp for
            for (std::int64_t i = 0; i < static_cast<std::int64_t>(num_cells()); ++i) {
                auto curpos = static_cast<pos_t<Dim>>(upos(i));
                auto pcell = cell(i);
                if (std::abs(pcell->crystallinity - 1.0) <= epsilon * (1.0 + pcell->crystallinity))
                    continue;
                
                Real delta = 0.0;
                std::map<grain_type*, grow_dir_type> nbhpcrysts_accdps;
                std::int64_t accdpmagn2 = 0;
                std::size_t numcrystednb = 0;
                for (auto nboff : get_nbhood_offset(i)) {
                    auto pnb = cell(nboff);
                    
                    if (pnb->crystallinity < 1.0 - epsilon * (1.0 + pcell->crystallinity) ||
                        pnb->grains.size() != 1) {
                        continue;
                    }

                    if (std::find(pcell->grains.begin(),
                                  pcell->grains.end(),
                                  pnb->grains.front()) == pcell->grains.end())
                        pcell->grains.push_back(pnb->grains.front());
                    
                    pos_t<Dim> deltapos = curpos - static_cast<pos_t<Dim>>(upos(nboff));
                    accdpmagn2 += deltapos.magnitude2();

                    nbhpcrysts_accdps[pnb->grains.front()] += static_cast<grow_dir_type>(deltapos);
                    ++numcrystednb;
                }

                for (auto& [pcryst, accdp] : nbhpcrysts_accdps) {
                    if (pcryst->material()->matproperty() == material_property::anisotropic) {
                        Real growth_factor = accdp.magnitude() / accdpmagn2;
                        for (auto& matergd : pcryst->material()->grow_dirs()) {
                            auto oriengd = spt::dot(pcryst->orientation().transposed(), matergd);
                            
                            delta += std::abs(spt::dot(accdp, oriengd)) * growth_factor;
                            // this will lead to alternative growth
                            //auto absdot = std::abs(spt::dot(accdp, oriengd)) * growth_factor;
                            //if (absdot > delta)
                            //    delta = absdot;
                        }
                    } else {
                        delta += accdp.magnitude() * 1.5;
                    }
                }
                Real factor = static_cast<Real>(std::clamp<std::int64_t>(
                    2 * accdpmagn2 - accdefdpmagn2, 0, accdefdpmagn2)) / accdefdpmagn2;
                delta += numcrystednb * factor;
                
                if (delta > epsilon)
                    m_cells_delta[i] += delta / default_nbhood_size() / 5;
            }

            #pragma omp barrier
            #pragma omp for
            for (std::int64_t i = 0; i < static_cast<std::int64_t>(num_cells()); ++i) {
                auto pcell = cell(i);
                pcell->crystallinity += m_cells_delta[i];
                m_cells_delta[i] = 0.0;

                if (pcell->crystallinity > 1.0 + epsilon * (1.0 + pcell->crystallinity)) {
                    pcell->crystallinity = 1.0;
                    update_nbhood_offsets_after_crystallization(i);
                }
            }
        }

        return true;
    }

    automata(std::size_t dimlen, shifts_fn<Dim> default_shift_fn, std::size_t default_range = 1)
        : automata(upos_t<Dim>::filled_with(dimlen), default_shift_fn, default_range) {}
    automata(const upos_t<Dim>& dimlens, shifts_fn<Dim> default_shift_fn, std::size_t default_range = 1) {
        m_default_shift_fn = default_shift_fn;
        set_dim_lens(dimlens);
        set_default_range(default_range);
        std::size_t new_num_cells = std::accumulate(m_dim_lens.x.begin(), m_dim_lens.x.end(),
                                                    static_cast<std::size_t>(1),
                                                    std::multiplies<std::size_t>());
        reserve(new_num_cells);
        m_nbhood_offsets.assign(new_num_cells, nbhood_offset_type());
        m_cells.assign(new_num_cells, nullptr);
        m_cells_delta.assign(num_cells(), 0.0);
    }


private:
    shifts_fn<Dim> m_default_shift_fn;
    std::size_t m_default_range;
    std::size_t m_default_nbhood_size;
    upos_t<Dim> m_dim_lens;

    cells_container m_cells;
    nbhood_offsets_container m_nbhood_offsets;

    cells_delta_container m_cells_delta;
    // todo: store common cells, such as empty cell and crystallized cell for each grain

    void set_default_nbhood_size() {
        m_default_nbhood_size = cgr::make_offsets<Dim>(m_default_shift_fn,
            0, m_dim_lens, m_default_range, std::nullopt).size();
    }
    void set_dim_lens(const upos_t<Dim>& dimlens) {
        m_dim_lens = dimlens;
    }

    bool is_nbhood_offsets_initialized = false;
    void initialize_nbhood_offsets() {
        #pragma omp parallel for
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(num_cells()); ++i)
            if (std::abs(cell(i)->crystallinity - 1.0)
                <= epsilon * (1.0 + cell(i)->crystallinity)) {
                for (auto nboff : make_offsets(i)) {
                    if (cell(nboff)->crystallinity <= epsilon &&
                        get_nbhood_offset(nboff).empty())
                        set_nbhood_offset(nboff);
                }
            }

        is_nbhood_offsets_initialized = true;
    }
    void update_nbhood_offsets_after_crystallization(std::size_t offset) {
        for (auto nboff : get_nbhood_offset(offset)) {
            if (cell(nboff)->crystallinity <= epsilon &&
                get_nbhood_offset(nboff).empty())
                set_nbhood_offset(nboff);
        }

        if (!get_nbhood_offset(offset).empty())
            clear_nbhood_offset(offset);
    }

    nbhood_offset_type make_offsets() const {
        upos_t<Dim> dimlens = upos_t<Dim>::filled_with(2 * m_default_range + 1);
        return cgr::make_offsets<Dim>(m_default_shift_fn,
            cgr::offset(pos_t<Dim>::filled_with(m_default_range), dimlens),
            dimlens, m_default_range, std::nullopt);
    }
};

} // namespace cgr
