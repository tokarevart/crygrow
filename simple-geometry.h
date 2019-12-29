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

    bool grains_includes_unsorted(std::vector<grain_type*> first, std::vector<grain_type*> second) const {
        std::sort(first.begin(), first.end());
        std::sort(second.begin(), second.end());
        return std::includes(first.begin(), first.end(), second.begin(), second.end());
    }
    bool grains_equals_unsorted(const std::vector<grain_type*>& first, const std::vector<grain_type*>& second) const {
        if (first.size() != second.size())
            return false;
        for (auto f : first)
            if (std::find(second.begin(), second.end(), f) == second.end())
                return false;
        return true;
    }

    offsets_container grain_boundaries(const cells_container& cells) const {
        offsets_container res;
        for (std::size_t i = 0; i < cells.size(); ++i)
            if (cells[i]->grains.size() > 1)
                res.push_back(i);
        return res;
    }
    std::vector<offsets_container> group_equal_grains_num_by_grains(const offsets_container& eqgrainsnum) const {
        std::vector<offsets_container> res;

    }
    std::vector<offsets_container> group_boundaries_by_grains_num(const offsets_container& bndoffsets) const {
        std::vector<offsets_container> res;
        for (auto off : bndoffsets) {
            std::size_t num_grains = m_automata->get_cell(off)->grains.size();
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
};

} // namespace cgr
