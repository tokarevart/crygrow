#pragma once
#include "automata-base.h"
#include <array>
#include <numeric>
#include <memory>


namespace cgr {

template <std::size_t Dim, typename Cell>
class vector_automata_base : automata_base<Dim, Cell> {
public:
    Cell* get_cell(const spt::vec<Dim, std::int64_t>& pos) const override final {
        auto actualpos = actual_pos(pos);
        for (std::size_t i = 0; i < Dim; i++)
            if (actualpos[i] >= m_dims_lens[i] ||
                actualpos[i] < 0)
                return nullptr;

        return m_cells[offset(actualpos)].get();
    }
    void reset_cell(const spt::vec<Dim, std::int64_t>& pos, Cell* ptr = nullptr) override final {
        auto actualpos = actual_pos(pos);
        bool need_resize = false;
        auto new_origin = m_origin;
        auto far_corner = m_origin + m_dims_lens;
        for (std::size_t i = 0; i < Dim; i++) {
            if (actualpos[i] >= m_dims_lens[i]) {
                far_corner[i] += 1 + actualpos[i] - m_dims_lens[i];
                far_corner[i] += (far_corner[i] - new_origin[i]) * 2;
                need_resize = true;
            }
            else if (actualpos[i] < 0) {
                new_origin[i] += actualpos[i];
                new_origin[i] += (new_origin[i] - far_corner[i]) * 2;
                need_resize = true;
            }
        }
        if (need_resize)
            resize_cells(new_origin, far_corner);

        m_cells[offset(actual_pos(pos))].reset(ptr);
    }

    virtual void resize_cells(const spt::vec<Dim, std::int64_t>& cuboid_corner0,
                              const spt::vec<Dim, std::int64_t>& cuboid_corner1) final {
        auto new_origin = cuboid_corner0;
        auto far_corner = cuboid_corner1;

        for (std::size_t i = 0; i < Dim; i++)
            if (new_origin[i] > far_corner[i])
                std::swap(new_origin[i], far_corner[i]);

        spt::vec<Dim, std::size_t> new_dims_lens = far_corner - new_origin;
        spt::vec<Dim, std::int64_t> dorigin = m_origin - new_origin;

        std::vector<std::unique_ptr<Cell>> new_cells(
            std::accumulate(new_dims_lens.x.begin(), new_dims_lens.x.end(), 1,
                            std::multiplies<std::size_t>()));

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

            new_i = offset(new_actualpos, new_dims_lens);
            new_cells[new_i] = std::move(m_cells[i]);
        }

        m_cells = std::move(new_cells);
        m_origin = new_origin;
        m_dims_lens = new_dims_lens;
    }
    virtual void reserve_cells(const spt::vec<Dim, std::int64_t>& cuboid_corner0,
                               const spt::vec<Dim, std::int64_t>& cuboid_corner1) final {
        auto new_origin = corner0;
        auto far_corner = corner1;

        for (std::size_t i = 0; i < Dim; i++)
            if (new_origin[i] > far_corner[i])
                std::swap(new_origin[i], far_corner[i]);

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

        resize_cells(new_origin, far_corner);

        //spt::vec<Dim, std::size_t> new_dims_lens = far_corner - new_origin;
        //spt::vec<Dim, std::size_t> dorigin = m_origin - new_origin;
        //
        //std::vector<std::unique_ptr<Cell>> new_cells(
        //    std::accumulate(new_dims_lens.x.begin(), new_dims_lens.x.end(), 1,
        //                    std::multiplies<std::size_t>()));
        //
        //for (std::size_t i = 0; i < m_cells.size(); i++) {
        //    auto new_actualpos = actual_pos(i) + dorigin;
        //    new_i = offset(new_actualpos, new_dims_lens);
        //    new_cells[new_i] = std::move(m_cells[i]);
        //}
        //
        //m_cells = std::move(new_cells);
        //m_origin = new_origin;
        //m_dims_lens = new_dims_lens;
    }
    virtual void shrink_to_fit_cells() final {
        // implement
    }

    virtual ~vector_automata_base() {}


private:
    spt::vec<Dim, std::int64_t> m_origin;
    spt::vec<Dim, std::size_t> m_dims_lens;
    std::vector<std::unique_ptr<Cell>> m_cells;

    std::size_t offset(const spt::vec<Dim, std::size_t>& pos) const {
        return offset(m_dims_lens);
    }
    std::size_t offset(const spt::vec<Dim, std::size_t>& pos, 
                       const spt::vec<Dim, std::size_t>& dims_lens) const {
        std::size_t res = pos.x[0];
        std::size_t mul = dims_lens[0];
        for (std::size_t i = 1; i < Dim; i++) {
            res += pos.x[i] * mul;
            mul *= dims_lens[i];
        }
        return res;
    }
    spt::vec<Dim, std::int64_t> actual_pos(const spt::vec<Dim, std::int64_t>& pos) const {
        return actual_pos(pos, m_origin);
    }
    spt::vec<Dim, std::int64_t> actual_pos(const spt::vec<Dim, std::int64_t>& pos,
                                         const spt::vec<Dim, std::int64_t>& origin) const {
        return pos - origin;
    }
    spt::vec<Dim, std::size_t> actual_pos(std::size_t offset) const {
        return actual_pos(offset, m_dims_lens);
    }
    spt::vec<Dim, std::size_t> actual_pos(std::size_t offset,
                                        const spt::vec<Dim, std::size_t>& dims_lens) const {
        // implement
    }
};

namespace gridimpl {
template <std::size_t Dim, typename Cell>
using vector = vector_automata_base<Dim, Cell>;
}

} // namespace cgr
