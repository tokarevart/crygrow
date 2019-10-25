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


template <std::size_t Dim, typename Cell, cell_mut_group CellMutGr>
class automata_base {
public:
    static constexpr std::size_t dim = Dim;
    using cell_type = Cell;
    using veci = spt::veci<Dim>;
    using vecu = spt::vecu<Dim>;
    using nbhood_type = nbhood<Dim, Cell>;
    // todo: using nbhood_offset_type = nbhood_offset<Dim>;
    using nbhood_pos_type = nbhood_pos<Dim>;
    using cells_container = std::vector<Cell*>;
    using nbhoods_pos_container = std::vector<nbhood_pos_type>;
    static constexpr cell_mut_group cell_mut_group = CellMutGr;
    
    std::size_t num_cells() const {
        return m_cells.size();
    }
    
    Cell* get_cell(std::size_t pos_offset) const {
        return m_cells[pos_offset];
    }
    Cell* get_cell(const veci& pos) const {
        return get_cell(offset(pos));
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
    
    Cell* try_get_cell(const veci& pos) const {
        return inside(pos) ? get_cell(offset(pos)) : nullptr;
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


    bool inside(const veci& pos) const {
        return cgr::inside(static_cast<vecu>(pos - m_origin), m_dim_lens);
    }
    veci pos(std::size_t pos_offset) const {
        return m_origin + upos(pos_offset);
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
        return cgr::offset(pos, m_dim_lens);
    }

    const nbhood_pos_type& get_nbhood_pos(std::size_t pos_offset) const {
        return m_nbhood_poses[pos_offset];
    }
    const nbhood_pos_type& get_nbhood_pos(const veci& pos) const {
        return inside(pos) ? get_nbhood_pos(offset(upos(pos))) : nullptr;
    }
    void  set_nbhood_pos(std::size_t pos_offset) {
        veci pos = this->pos(pos_offset);
        if (inside(pos)) {
            m_nbhood_poses[pos_offset]
                = make_nbhood_pos<Dim>(pos, default_nbhood_kind(), default_range(),
                                       [this](const veci& pos) -> bool { return this->try_get_cell(pos); });
        }
    }
    void  set_nbhood_pos(const veci& pos) {
        if (inside(pos)) {
            m_nbhood_poses[offset(pos)] 
                = make_nbhood_pos<Dim>(pos, default_nbhood_kind(), default_range(),
                                       [this](const veci& pos) -> bool { return try_get_cell(pos); });
        }
    }
    void  set_nbhood_pos(const veci& pos, nbhood_pos_type&& nbhpos) {
        if (inside(pos))
            m_nbhood_poses[offset(pos)] = std::move(nbhpos);
    }
    nbhood_type get_nbhood(const veci& pos) const {
        return inside(pos) ? make_nbhood<Dim, Cell>(
            get_nbhood_pos(offset(pos)),
            [this](const veci& pos) -> Cell* { return try_get_cell(pos); }) : nbhood_type();
    }
    
    void reserve(std::size_t count) {
        reserve_cells(count);
        reserve_nbhoods_poses(count);
    }
    void reserve_cells(std::size_t count) {
        m_cells.reserve(count);
    }
    void reserve_nbhoods_poses(std::size_t count) {
        m_nbhood_poses.reserve(count);
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
    nbhood_kind default_nbhood_kind() const {
        return m_default_nbhood_kind;
    }
    std::size_t default_nbhood_size() const {
        return m_default_nbhood_size;
    }

    void set_defaults(std::size_t range, nbhood_kind kind) {
        m_default_range = range;
        m_default_nbhood_kind = kind;
        set_default_nbhood_size();
    }
    void set_default_range(std::size_t range) {
        m_default_range = range;
        set_default_nbhood_size();
    }
    void set_default_nbhood_kind(nbhood_kind kind) {
        m_default_nbhood_kind = kind;
        set_default_nbhood_size();
    }

    virtual bool stop_condition() const = 0;
    virtual bool iterate() = 0;

    automata_base(const veci& corner0, const veci& corner1, 
                  std::size_t default_range, nbhood_kind default_nbhood_kind) {
        set_defaults(default_range, default_nbhood_kind);
        set_origin_and_dim_lens(corner0, corner1);
        std::size_t new_num_cells = std::accumulate(m_dim_lens.x.begin(), m_dim_lens.x.end(), 
                                                    static_cast<std::size_t>(1),
                                                    std::multiplies<std::size_t>());
        reserve(new_num_cells);
        m_nbhood_poses.assign(new_num_cells, nbhood_pos_type());
        m_cells.assign(new_num_cells, nullptr);
    }

    virtual ~automata_base() {}


private:
    std::size_t m_default_range;
    nbhood_kind m_default_nbhood_kind;
    std::size_t m_default_nbhood_size;

    veci m_origin;
    vecu m_dim_lens;

    cells_container m_cells;
    nbhoods_pos_container m_nbhood_poses;

    void set_default_nbhood_size() {
        m_default_nbhood_size = make_nbhood_pos<Dim>({ 0, 0 }, default_nbhood_kind(), default_range()).size();
    }
    
    bool inside(const vecu& pos) const {
        return inside(pos, m_dim_lens);
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
