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
#include "simple-cell.h"


namespace cgr {

template <std::size_t Dim, nbhood_kind NbhoodKind = nbhood_kind::euclid, typename Real = double>
class simple_automata {
public:
    static constexpr std::size_t dim = Dim;
    static constexpr cgr::nbhood_kind nbhood_kind = NbhoodKind;
    static constexpr Real epsilon = std::numeric_limits<Real>::epsilon();
    using veci = spt::veci<Dim>;
    using vecu = spt::vecu<Dim>;
    using cell_type = cgr::simple_cell<Dim, Real>;
    using cells_container = std::vector<cell_type*>;
    using nbhood_type = cgr::nbhood<Dim, cell_type>;
    using nbhood_pos_type = cgr::nbhood_pos<Dim>;
    using nbhood_offset_type = cgr::nbhood_offset;
    using nbhood_offsets_container = std::vector<nbhood_offset_type>;
    using cells_delta_container = std::vector<Real>;
    using grain_type = typename cell_type::grain_type;
    using material_type = typename grain_type::material_type;
    using orientation_type = typename grain_type::orientation_type;
    using grow_dir = typename material_type::grow_dir;

    std::size_t num_cells() const {
        return m_cells.size();
    }
    const cells_container& cells() const {
        return m_cells;
    }

    cell_type* get_cell(std::size_t offset) const {
        return m_cells[offset];
    }
    cell_type* get_cell(const vecu& pos) const {
        return get_cell(offset(pos));
    }
    cell_type* get_cell(const veci& pos) const {
        return get_cell(offset(pos));
    }
    void  set_cell(std::size_t offset, const cell_type* new_cell) {
        m_cells[offset] = const_cast<cell_type*>(new_cell);
    }
    void  set_cell(const veci& pos, const cell_type* new_cell) {
        m_cells[offset(pos)] = const_cast<cell_type*>(new_cell);
    }
    template <typename PosesContainer, typename CellsContainer>
    void  set_cells(const PosesContainer& poses, const CellsContainer& cells) {
        auto pos_it = poses.begin();
        auto cell_it = cells.begin();
        for (; pos_it != poses.end(); ++pos_it, ++cell_it)
            set_cell(*pos_it, &(*cell_it));
    }

    cell_type* try_get_cell(const vecu& pos) const {
        return inside(pos) ? get_cell(pos) : nullptr;
    }
    cell_type* try_get_cell(const veci& pos) const {
        return inside(pos) ? get_cell(pos) : nullptr;
    }
    bool  try_set_cell(const veci& pos, const cell_type* new_cell) {
        if (inside(pos)) {
            m_cells[offset(pos)] = const_cast<cell_type*>(new_cell);
            return true;
        }

        return false;
    }
    template <typename PosesContainer, typename CellsContainer>
    void  try_set_cells(const PosesContainer& poses, const CellsContainer& cells) {
        auto pos_it = poses.begin();
        auto cell_it = cells.begin();
        for (; pos_it != poses.end(); ++pos_it, ++cell_it)
            try_set_cell(*pos_it, &(*cell_it));
    }

    bool inside(const vecu& pos) const {
        return inside(pos, m_dim_lens);
    }
    bool inside(const veci& pos) const {
        for (auto& e : pos.x)
            if (e < 0)
                return false;
        return cgr::inside(static_cast<vecu>(pos), m_dim_lens);
    }
    vecu upos(std::size_t offset) const {
        return cgr::upos(offset, m_dim_lens);
    }
    std::size_t offset(const vecu& pos) const {
        return cgr::offset(static_cast<veci>(pos), m_dim_lens);
    }

    nbhood_offset_type make_nbhood_offset() const {
        return cgr::make_nbhood_offset<NbhoodKind, Dim>(
            cgr::offset(veci::filled_with(m_default_range), vecu::filled_with(2 * m_default_range + 1)),
            vecu::filled_with(2 * m_default_range + 1), m_default_range);
    }
    nbhood_offset_type make_nbhood_offset(std::size_t pos_offset) const {
        return cgr::make_nbhood_offset<NbhoodKind, Dim>(
            pos_offset, m_dim_lens, m_default_range,
            [this](const veci& pos) -> bool { return try_get_cell(pos); });
    }
    const nbhood_offset_type& get_nbhood_offset(std::size_t pos_offset) const {
        return m_nbhood_offsets[pos_offset];
    }
    void set_nbhood_offset(std::size_t pos_offset) {
        m_nbhood_offsets[pos_offset] = make_nbhood_offset(pos_offset);
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

    void erase(const veci& pos) {
        erase_cell(pos);
        erase_nbhood_poses(pos);
    }
    void erase_cell(const veci& pos) {
        set_cell(pos, nullptr);
    }
    void erase_nbhood_poses(const veci& pos) {
        if (inside(pos))
            m_nbhood_poses[offset(pos)].clear();
    }

    vecu dim_lens() const {
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
            if (get_cell(i)->crystallinity < 1.0 - epsilon * (1.0 + get_cell(i)->crystallinity))
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

        auto defnbhoffset = make_nbhood_offset();
        std::int64_t accdefdpmagn2 = 0;
        for (auto nboff : defnbhoffset)
            accdefdpmagn2 += (
                veci::filled_with(default_range())
                - cgr::upos(nboff, vecu::filled_with(2 * default_range() + 1))
                ).magnitude2();
        
        #pragma omp parallel 
        {
            #pragma omp for
            for (std::int64_t i = 0; i < static_cast<std::int64_t>(num_cells()); ++i) {
                auto curpos = static_cast<veci>(upos(i));
                auto pcell = get_cell(i);
                if (std::abs(pcell->crystallinity - 1.0) <= epsilon * (1.0 + pcell->crystallinity))
                    continue;
                
                Real delta = 0.0;
                std::map<grain_type*, grow_dir> nbhpcrysts_accdps;
                std::int64_t accdpmagn2 = 0;
                std::size_t numcrystednb = 0;
                for (auto nboff : get_nbhood_offset(i)) {
                    auto pnb = get_cell(nboff);
                    
                    if (pnb->crystallinity < 1.0 - epsilon * (1.0 + pcell->crystallinity) ||
                        pnb->grains.size() != 1) {
                        continue;
                    }

                    if (std::find(pcell->grains.begin(),
                                  pcell->grains.end(),
                                  pnb->grains.front()) == pcell->grains.end())
                        pcell->grains.push_back(pnb->grains.front());
                    
                    veci deltapos = curpos - static_cast<veci>(upos(nboff));
                    accdpmagn2 += deltapos.magnitude2();

                    nbhpcrysts_accdps[pnb->grains.front()] += static_cast<grow_dir>(deltapos);
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
                        delta += accdp.magnitude();
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
                auto pcell = get_cell(i);
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

    simple_automata(std::size_t dimlen, std::size_t default_range = 1)
        : simple_automata(vecu::filled_with(dimlen), default_range) {}
    simple_automata(const vecu& dimlens, std::size_t default_range = 1) {
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
    std::size_t m_default_range;
    std::size_t m_default_nbhood_size;
    vecu m_dim_lens;

    cells_container m_cells;
    nbhood_offsets_container m_nbhood_offsets;

    cells_delta_container m_cells_delta;
    // todo: store common cells, such as empty cell and crystallized cell for each crystallite

    void set_default_nbhood_size() {
        m_default_nbhood_size = cgr::make_nbhood_offset<NbhoodKind, Dim>(
            0, m_dim_lens, m_default_range).size();
    }
    void set_dim_lens(const vecu& dimlens) {
        m_dim_lens = dimlens;
    }

    bool is_nbhood_offsets_initialized = false;
    void initialize_nbhood_offsets() {
        #pragma omp parallel for
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(num_cells()); ++i)
            if (std::abs(get_cell(i)->crystallinity - 1.0)
                <= epsilon * (1.0 + get_cell(i)->crystallinity)) {
                for (auto nboff : make_nbhood_offset(i)) {
                    if (get_cell(nboff)->crystallinity <= epsilon &&
                        get_nbhood_offset(nboff).empty())
                        set_nbhood_offset(nboff);
                }
            }

        is_nbhood_offsets_initialized = true;
    }
    void update_nbhood_offsets_after_crystallization(std::size_t offset) {
        for (auto nboff : get_nbhood_offset(offset)) {
            if (get_cell(nboff)->crystallinity <= epsilon &&
                get_nbhood_offset(nboff).empty())
                set_nbhood_offset(nboff);
        }

        if (!get_nbhood_offset(offset).empty())
            clear_nbhood_offset(offset);
    }
};

} // namespace cgr
