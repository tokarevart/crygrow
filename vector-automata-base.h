// Copyright © 2019 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <array>
#include <vector>
#include <numeric>
#include <memory>
#include "pos-automata-base.h"
#include "cell-mutability.h"


namespace cgr {

template <std::size_t Dim, typename Cell, cell_mut_group CellMutGr>
class vector_automata_base : pos_automata_base<Dim, Cell> {
public:
    using cells_container = std::vector<std::unique_ptr<Cell>>;
    using iterator = cell_iterator<vector_automata_base<Dim, Cell, CellMutGr>>;
    static constexpr cell_mut_group cell_mut_group = CellMutGr;

    iterator begin() const {
        return { m_cells.cbegin() };
    }
    iterator end() const {
        return { m_cells.cend() };
    }

    Cell* cell(const veci& pos) const override {
        return inside(pos) ? m_cells[offset(actual_pos(pos))].get() : nullptr;
    }
    void  cell(const veci& pos, const Cell* new_cell) override {
        if (!new_cell) {
            if (auto pcell = cell(pos); pcell)
                pos_automata_base<Dim, Cell>::erase(pcell);
            if (inside(pos))
                m_cells[offset(actual_pos(pos))].reset(new_cell);
        } else {
            try_reallocate(pos);
            m_cells[offset(actual_pos(pos))].reset(new_cell);
        }
    }

    // NOTE: don't recomend reduce the sizes along the axes for now
    // NOTE: needs implemented actual_pos(offset, dims_lens)
    void resize(const veci& corner0, const veci& corner1) {
        auto new_origin = corner0;
        auto far_corner = corner1;

        for (std::size_t i = 0; i < Dim; i++)
            sort2(new_origin[i], far_corner[i]);
        
        vecu new_dims_lens = far_corner + veci{ 1, 1, 1 } - new_origin;
        veci dorigin = m_origin - new_origin;

        auto num_new_cells = std::accumulate(new_dims_lens.x.begin(), new_dims_lens.x.end(), 1,
                                             std::multiplies<std::size_t>());
        pos_automata_base<Dim, Cell>::reserve(num_new_cells);
        std::vector<std::unique_ptr<Cell>> new_cells(num_new_cells);

        for (std::size_t i = 0; i < m_cells.size(); i++) {
            auto new_actualpos = actual_pos(i) + dorigin;
            bool abroad = false;
            for (std::size_t i = 0; i < Dim; i++)
                if (new_actualpos[i] < 0 || 
                    new_actualpos[i] >= new_dims_lens[i]) {
                    abroad = true;
                    break;
                }
            if (abroad)
                continue;

            auto new_i = offset(new_actualpos, new_dims_lens);
            new_cells[new_i] = std::move(m_cells[i]);
        }

        m_cells = std::move(new_cells);
        m_origin = new_origin;
        m_dims_lens = new_dims_lens;
    }
    void reserve(const veci& corner0, const veci& corner1) {
        auto new_origin = corner0;
        auto far_corner = corner1;

        for (std::size_t i = 0; i < Dim; i++)
            sort2(new_origin[i], far_corner[i]);

        bool need_resize = false;
        for (std::size_t i = 0; i < Dim; i++)
            if (new_origin[i] < m_origin[i] ||
                far_corner[i] > m_origin[i] + m_dims_lens[i]) {
                need_resize = true;
                break;
            }
        if (!need_resize)
            return;

        for (std::size_t i = 0; i < Dim; i++)
            new_origin[i] = std::min(new_origin[i], m_origin[i]);

        for (std::size_t i = 0; i < Dim; i++)
            far_corner[i] = std::max(far_corner[i], m_origin[i] + m_dims_lens[i]);

        resize(new_origin, far_corner);
    }
    void shrink_to_fit() {
        // implementation
    }

    virtual ~vector_automata_base() {}


private:
    veci m_origin;
    vecu m_dims_lens;
    cells_container m_cells;

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
    bool try_reallocate(const veci& pos) {
        auto actualpos = actual_pos(pos);
        bool need_resize = false;
        auto new_origin = m_origin;
        auto far_corner = m_origin + m_dims_lens;
        for (std::size_t i = 0; i < Dim; i++) {
            if (actualpos[i] >= m_dims_lens[i]) {
                far_corner[i] += 1 + actualpos[i] - m_dims_lens[i];
                far_corner[i] += (far_corner[i] - new_origin[i]) * 2;
                need_resize = true;
            } else if (actualpos[i] < 0) {
                new_origin[i] += actualpos[i];
                new_origin[i] += (new_origin[i] - far_corner[i]) * 2;
                need_resize = true;
            }
        }
        if (need_resize)
            resize(new_origin, far_corner);

        return need_resize;
    }
    void sort2(veci::value_type& first, veci::value_type& second) {
        if (first > second)
            std::swap(first, second);
    }
};


template <std::size_t Dim, typename Cell, cell_mut_group CellMutGr>
class cell_iterator<vector_automata_base<Dim, Cell, CellMutGr>>
    : cell_iterator_base<vector_automata_base<Dim, Cell, CellMutGr>> {
public:
    cell_type* to_ptr() const {
        return m_it->get();
    }
    cell_iterator next_until(cell_iterator end) {
        return cell_iterator(*this).adv_next_until(end);
    }
    cell_iterator& adv_next_until(cell_iterator end) {
        if (*this == end)
            return end;

        ++m_it;
        if constexpr (std::is_same_v<CellMutGr, cell_mut_group::universal>)
            return to_ptr()->mutability() == cell_mut::constant_cell ? adv_next_until(end) : *this;
        else if constexpr (std::is_same_v<CellMutGr, cell_mut_group::mutable_only>)
            return *this;
        else
            static_assert(false);
    }

    bool operator==(const cell_iterator& other) const {
        return m_it == other.m_it;
    }
    bool operator!=(const cell_iterator& other) const {
        return m_it != other.m_it;
    }
    cell_type& operator*() const {
        return *to_ptr();
    }

    cell_iterator(from_iterator it) : m_it{it} {}


private:
    from_iterator m_it;
};

} // namespace cgr
