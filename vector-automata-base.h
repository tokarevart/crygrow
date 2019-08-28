#pragma once
#include "automata-base.h"
#include <array>
#include <numeric>
#include <memory>


namespace cgr {

template <std::size_t Dim, typename Cell>
class vector_automata_base 
    : automata_base<Dim, Cell> {
public:
    using cells_container_type = std::vector<std::unique_ptr<Cell>>;
    using iterator = cell_iterator<vector_automata_base<Dim, Cell>>;

    iterator begin() const {
        return { m_cells.begin() };
    }
    iterator end() const {
        return { m_cells.end() };
    }

    Cell* get(const veci& pos) const override final {
        auto actualpos = actual_pos(pos);
        for (std::size_t i = 0; i < Dim; i++)
            if (actualpos[i] >= m_dims_lens[i] ||
                actualpos[i] < 0)
                return nullptr;

        return m_cells[offset(actualpos)].get();
    }
    void reset(const veci& pos, Cell* ptr = nullptr) override final {
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
            resize(new_origin, far_corner);

        m_cells[offset(actual_pos(pos))].reset(ptr);
    }

    virtual void resize(const veci& corner0,
                        const veci& corner1) final {
        auto new_origin = corner0;
        auto far_corner = corner1;

        for (std::size_t i = 0; i < Dim; i++)
            sort2(new_origin[i], far_corner[i]);

        vecu new_dims_lens = far_corner - new_origin;
        veci dorigin = m_origin - new_origin;

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
    virtual void reserve(const veci& corner0,
                         const veci& corner1) final {
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
    virtual void shrink_to_fit() final {
        // implementation
    }

    virtual ~vector_automata_base() {}


private:
    veci m_origin;
    vecu m_dims_lens;
    cells_container_type m_cells;

    std::size_t offset(const vecu& pos) const {
        return offset(pos, m_dims_lens);
    }
    std::size_t offset(const vecu& pos,
                       const vecu& dims_lens) const {
        std::size_t res = pos.x[0];
        std::size_t mul = dims_lens[0];
        for (std::size_t i = 1; i < Dim; i++) {
            res += pos.x[i] * mul;
            mul *= dims_lens[i];
        }
        return res;
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
    template <typename T>
    void sort2(T& first, T& second) {
        if (first > second)
            std::swap(first, second);
    }
};


template <std::size_t Dim, typename Cell>
class cell_iterator<vector_automata_base<Dim, Cell>> : cell_iterator_base<vector_automata_base<Dim, Cell>> {
public:
    // implementation

    cell_iterator(from_iterator it) 
        : m_it{ it } {}


private:
    from_iterator m_it;
};

} // namespace cgr
