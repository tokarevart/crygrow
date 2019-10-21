// Copyright © 2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <memory>
#include <unordered_map>
#include <optional>
#include "neighborhood.h"
#include "cell-mutability.h"
#include "vec.h"


namespace cgr {

template <typename Automata>
class cell_iterator;


template <std::size_t Dim, typename Cell, cell_mut_group CellMutGr>
class automata_base {
public:
    static constexpr std::size_t dim = Dim;
    using cell_type = Cell;
    using veci = spt::vec<Dim, std::int64_t>;
    using vecu = spt::vec<Dim, std::uint64_t>;
    using nbhood_type = nbhood<Dim, Cell>;
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
        return inside(pos) ? get_cell(offset(actual_pos(pos))) : nullptr;
    }
    void  set_cell(const veci& pos, const Cell* new_cell) {
        if (inside(pos))
            m_cells[offset(actual_pos(pos))] = new_cell;
    }
    void  set_cells(const std::vector<veci>& poses, const std::vector<Cell*>& cells) {
        for (std::size_t i = 0; i < poses.size(); i++)
            set_cell(poses[i], cells[i]);
    }
    void  set_cells(const std::vector<veci>& poses, const std::vector<Cell>& cells) {
        for (std::size_t i = 0; i < poses.size(); i++)
            set_cell(poses[i], &cells[i]);
    }

    veci get_pos(std::size_t pos_offset) const {
        return origin() + actual_pos(i);
    }

    const nbhood_pos_type& get_nbhood_pos(std::size_t pos_offset) const {
        return m_nbhoods_poses[pos_offset];
    }
    const nbhood_pos_type& get_nbhood_pos(const veci& pos) const {
        return inside(pos) ? get_nbhood_pos(offset(actual_pos(pos))) : nullptr;
    }
    void  set_nbhood_pos(const veci& pos) {
        if (inside(pos)) {
            m_nbhoods_poses[offset(actual_pos(pos))] 
                = make_nbhood_pos(pos, default_nbhood_kind(), default_range(),
                                  [this](const veci& pos) -> bool { return get_cell(pos); });
        }
    }
    void  set_nbhood_pos(const veci& pos, nbhood_pos_type&& nbhpos) {
        if (inside(pos))
            m_nbhoods_poses[offset(actual_pos(pos))] = std::move(nbhpos);
    }
    nbhood_type get_nbhood(const veci& pos) const {
        return inside(pos) ? make_nbhood(get_nbhood_pos(offset(actual_pos(pos))),
                                         [this](const veci& pos) { return get_cell(pos); }) : nullptr;
    }
    
    void reserve(std::size_t count) {
        reserve_cells(count);
        reserve_nbhoods_poses(count);
    }
    void reserve_cells(std::size_t count) {
        m_cells.reserve(count);
    }
    void reserve_nbhoods_poses(std::size_t count) {
        m_nbhoods_poses.reserve(count);
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
            m_nbhoods_poses[offset(actual_pos(pos))].clear();
    }
    
    std::size_t default_range() const {
        return m_default_range;
    }
    nbhood_kind default_nbhood_kind() const {
        return m_default_nbhood_kind;
    }
    
    virtual bool stop_condition() const = 0;
    virtual bool iterate() = 0;

    automata_base(const veci& corner0, const veci& corner1, 
                  std::size_t default_range, nbhood_kind default_nbhood_kind)
        : m_default_range{default_range}, m_default_nbhood_kind{default_nbhood_kind} {
        veci new_origin = corner0;
        veci far_corner = corner1;

        for (std::size_t i = 0; i < Dim; i++)
            sort2(new_origin[i], far_corner[i]);
        m_origin = new_origin;

        m_dims_lens = far_corner + veci{ 1, 1, 1 } - new_origin;
        std::size_t new_num_cells = std::accumulate(m_dims_lens.x.begin(), m_dims_lens.x.end(), 1,
                                                    std::multiplies<std::size_t>());
        reserve(new_num_cells);
        m_nbhoods_poses.assign(new_num_cells, nbhood_pos_type);
        m_cells.assign(new_num_cells, nullptr);
    }

    virtual ~automata_base() {}


private:
    nbhood_kind m_default_nbhood_kind;
    std::size_t m_default_range;

    veci m_origin;
    vecu m_dims_lens;

    cells_container m_cells;
    nbhoods_pos_container m_nbhoods_poses;

    veci origin() const {
        return m_origin;
    }
    vecu dims_lens() const {
        return m_dims_lens;
    }

    std::size_t offset(const vecu& pos) const {
        return offset(pos, m_dims_lens);
    }
    std::size_t offset(const vecu& pos, const vecu& dims_lens) const {
        std::size_t res = pos.x[0];
        std::size_t mul = dims_lens[0];
        for (std::size_t i = 1; i < Dim; i++) {
            res += pos.x[i] * mul;
            mul *= dims_lens[i];
        }
        return res;
    }

    bool inside(const veci& pos) const {
        return inside(pos, m_dims_lens);
    }
    bool inside(const veci& pos, const vecu& dims_lens) const {
        return inside(pos, dims_lens, m_origin);
    }
    bool inside(const veci& pos, const vecu& dims_lens, const veci& origin) const {
        return inside(actual_pos(pos, origin), dims_lens);
    }
    bool inside(const vecu& pos) const {
        return inside(pos, m_dims_lens);
    }
    bool inside(const vecu& pos, const vecu& dims_lens) const {
        for (std::size_t i = 0; i < dims_lens.dim; i++)
            if (pos[i] >= dims_lens[i])
                return false;

        return true;
    }

    veci actual_pos(const veci& pos) const {
        return actual_pos(pos, m_origin);
    }
    veci actual_pos(const veci& pos, const veci& origin) const {
        return pos - origin;
    }
    vecu actual_pos(std::size_t offset) const {
        return actual_pos(offset, m_dims_lens);
    }
    vecu actual_pos(std::size_t offset, const vecu& dims_lens) const {
        // implementation
    }

    void sort2(veci::value_type& first, veci::value_type& second) {
        if (first > second)
            std::swap(first, second);
    }
};

} // namespace cgr
