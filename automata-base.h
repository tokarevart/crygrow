// Copyright © 2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <memory>
#include <optional>
#include "neighborhood.h"
#include "cell-mutability.h"
#include "vec.h"
#include "sptops.h"
#include "cgralgs.h"


namespace cgr {

template <typename Automata>
class cell_iterator;


template <std::size_t Dim, nbhood_kind NbhoodKind, typename Cell, cell_mut_group CellMutGr>
class automata_base {
public:
    static constexpr std::size_t dim = Dim;
    using cell_type = Cell;
    using veci = spt::veci<Dim>;
    using vecu = spt::vecu<Dim>;
    using nbhood_type = nbhood<Dim, Cell>;
    using nbhood_offset_type = nbhood_offset;
    using nbhood_pos_type = nbhood_pos<Dim>;
    using cells_container = std::vector<Cell*>;
    using nbhood_offsets_container = std::vector<nbhood_offset_type>;
    static constexpr cgr::nbhood_kind nbhood_kind = NbhoodKind;
    static constexpr cgr::cell_mut_group cell_mut_group = CellMutGr;

    std::size_t num_cells() const {
        return m_cells.size();
    }

    Cell* get_cell(std::size_t offset) const {
        return m_cells[offset];
    }
    Cell* get_cell(const vecu& pos) const {
        return get_cell(offset(pos));
    }
    Cell* get_cell(const veci& pos) const {
        return get_cell(offset(pos));
    }
    void  set_cell(std::size_t offset, const Cell* new_cell) {
        m_cells[offset] = const_cast<Cell*>(new_cell);
    }
    void  set_cell(const veci& pos, const Cell* new_cell) {
        m_cells[offset(pos)] = const_cast<Cell*>(new_cell);
    }
    template <typename PosesContainer, typename CellsContainer>
    void  set_cells(const PosesContainer& poses, const CellsContainer& cells) {
        auto pos_it = poses.begin();
        auto cell_it = cells.begin();
        for (; pos_it != poses.end(); ++pos_it, ++cell_it)
            set_cell(*pos_it, &(*cell_it));
    }

    Cell* try_get_cell(const vecu& pos) const {
        return inside(pos) ? get_cell(pos) : nullptr;
    }
    Cell* try_get_cell(const veci& pos) const {
        return inside(pos) ? get_cell(pos) : nullptr;
    }
    bool  try_set_cell(const veci& pos, const Cell* new_cell) {
        if (inside(pos)) {
            m_cells[offset(pos)] = const_cast<Cell*>(new_cell);
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
        return cgr::inside(static_cast<vecu>(pos - m_origin), m_dim_lens);
    }
    veci pos(std::size_t offset) const {
        return m_origin + upos(offset);
    }
    vecu upos(std::size_t offset) const {
        return cgr::upos(offset, m_dim_lens);
    }
    vecu upos(const veci& pos) const {
        return pos - m_origin;
    }
    std::size_t offset(const veci& pos) const {
        return offset(upos(pos));
    }
    std::size_t offset(const vecu& pos) const {
        return cgr::offset(static_cast<veci>(pos), m_dim_lens);
    }

    nbhood_offset_type make_nbhood_offset() const {
        return cgr::make_nbhood_offset<NbhoodKind, Dim>(
            0, m_dim_lens, m_default_range);
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

    virtual bool stop_condition() const = 0;
    virtual bool iterate() = 0;

    automata_base(const veci& corner0, const veci& corner1, 
                  std::size_t default_range) {
        set_origin_and_dim_lens(corner0, corner1);
        set_default_range(default_range);
        std::size_t new_num_cells = std::accumulate(m_dim_lens.x.begin(), m_dim_lens.x.end(), 
                                                    static_cast<std::size_t>(1),
                                                    std::multiplies<std::size_t>());
        reserve(new_num_cells);
        m_nbhood_offsets.assign(new_num_cells, nbhood_offset_type());
        m_cells.assign(new_num_cells, nullptr);
    }

    virtual ~automata_base() {}


private:
    std::size_t m_default_range;
    std::size_t m_default_nbhood_size;

    veci m_origin;
    vecu m_dim_lens;

    cells_container m_cells;
    nbhood_offsets_container m_nbhood_offsets;

    void set_default_nbhood_size() {
        m_default_nbhood_size = cgr::make_nbhood_offset<NbhoodKind, Dim>(
            0, m_dim_lens, m_default_range).size();
    }    
    void set_origin_and_dim_lens(const veci& corner0, const veci& corner1) {
        veci new_origin = corner0;
        veci far_corner = corner1;

        spt::sort_elementwise(new_origin, far_corner);
        m_origin = new_origin;
        m_dim_lens = far_corner - new_origin;
    }
};

} // namespace cgr
