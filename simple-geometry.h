#pragma once
#include "simple-automata.h"


namespace cgr {

template <nbhood_kind NbhoodKind = nbhood_kind::euclid, typename Real = double>
class simple_geometry_tool {
public:
    static constexpr std::size_t dim = 3;
    using automata_type = cgr::simple_automata<dim, NbhoodKind, Real>;
    using cell_type = typename automata_type::cell_type;
    using cells_container = typename automata_type::cells_container;
    using offsets_container = std::vector<std::size_t>;
    using grain_type = typename automata_type::grain_type;
    using grains_container = typename cell_type::grains_container;
    using grains_num_then_boundaries = std::vector<std::vector<grains_container>>;

    bool grains_includes_unsorted(const grains_container& first, const grains_container& second) const {
        if (first.size() < second.size())
            return false;
        for (auto s : second)
            if (std::find(first.begin(), first.end(), s) == first.end())
                return false;
        return true;
    }
    // each grain of both vectors is unique
    bool grains_equal_unsorted(const grains_container& first, const grains_container& second) const {
        if (first.size() != second.size())
            return false;
        for (auto f : first)
            if (std::find(second.begin(), second.end(), f) == second.end())
                return false;
        return true;
    }
    std::vector<offsets_container>::iterator grains_find_in_offsets(
        std::vector<offsets_container>& offsets, const grains_container& grains) const {
        auto it = offsets.begin();
        for (; it < offsets.end(); ++it)
            if (grains_equal_unsorted(get_grains(it->front()), grains))
                break;
        return it;
    }

    bool is_surface_joint(const grains_container& bndgrains, const grains_num_then_boundaries& grconts) const {
        return bndgrains.size() == 2;
    }
    bool is_linear_joint(const grains_container& bndgrains, const grains_num_then_boundaries& grconts) const {
        if (bndgrains.size() < 3)
            return false;
        //
    }
    bool is_point_joint(const grains_container& bndgrains, const grains_num_then_boundaries& grconts) const {
        if (bndgrains.size() < 4)
            return false;
        return !is_linear_joint(bndgrains, grconts);
    }

    offsets_container grain_boundaries() const {
        offsets_container res;
        for (std::size_t i = 0; i < num_cells(); ++i)
            if (get_grains(i).size() > 1)
                res.push_back(i);
        return res;
    }
    const grains_container& boundary_grains(const offsets_container& boundary) const {
        return get_grains(boundary.front());
    }
    std::vector<grains_container> boundary_grains(const std::vector<offsets_container>& boundaries) const {
        std::vector<std::vector<grain_type*>> res;
        res.reserve(boundaries.size());
        for (auto& offcont : boundaries)
            res.push_back(boundary_grains(offcont));
        return res;
    }
    std::vector<offsets_container> group_boundaries_by_grains(const offsets_container& bndoffsets) const {
        std::vector<offsets_container> res;
        for (auto off : bndoffsets) {
            auto found = grains_find_in_offsets(res, get_grains(off));
            if (found == res.end()) {
                res.push_back(offsets_container{ off });
            } else {
                found->push_back(off);
            }
        }
        return res;
    }
    std::vector<offsets_container> group_boundaries_by_grains_num(const offsets_container& bndoffsets) const {
        std::vector<offsets_container> res;
        for (auto off : bndoffsets) {
            std::size_t num_grains = get_grains(off).size();
            if (num_grains > res.size())
                while (res.size() < num_grains)
                    res.emplace_back();

            res[num_grains].push_back(off);
        }        
        return res;
    }

    simple_geometry_tool(const automata_type* automata)
        : m_automata(automata) {}


private:
    const automata_type* m_automata;

    std::size_t num_cells() const {
        return m_automata->num_cells();
    }
    cell_type* get_cell(std::size_t offset) const {
        return m_automata->get_cell(offset);
    }
    const grains_container& get_grains(std::size_t offset) const {
        return get_cell(offset)->grains;
    }
};

} // namespace cgr
