#pragma once
#include "simple-automata.h"


namespace cgr {

template <nbhood_kind NbhoodKind = nbhood_kind::euclid, typename Real = double>
class simple_geometry_tool {
public:
    static constexpr std::size_t dim = 3;
    using automata_type = cgr::simple_automata<dim, NbhoodKind, Real>;
    using pos_type = spt::vec3<Real>;
    using cell_type = typename automata_type::cell_type;
    using cells_container = typename automata_type::cells_container;
    using offsets_container = std::vector<std::size_t>;
    using grain_type = typename automata_type::grain_type;
    using grains_container = typename cell_type::grains_container;
    // grains num then boundaries sorted
    using grouped2_grains_container = std::vector<std::vector<grains_container>>;
    using grouped2_offsets_container = std::vector<std::vector<offsets_container>>;

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

    offsets_container offsets_intersection(const offsets_container& first, const offsets_container& second) const {
        //
    }

    // pjoint is a point joint (or point boundary)
    bool is_supreme_pjoint(const grains_container& bndgrains, 
                           const grouped2_grains_container& grconts) const {
        if (bndgrains.size() < 4)
            return false;

        for (auto it = grconts.begin() + bndgrains.size(); it < grconts.end(); ++it)
            for (auto& grcont : *it)
                if (grains_includes_unsorted(grcont, bndgrains))
                    return false;

        return true;
    }

    offsets_container boundaries_offsets() const {
        offsets_container res;
        for (std::size_t i = 0; i < num_cells(); ++i)
            if (get_grains(i).size() > 1)
                res.push_back(i);
        return res;
    }
    const grains_container& boundary_grains(const offsets_container& boundary) const {
        return get_grains(boundary.front());
    }
    std::vector<grains_container> boundaries_grains(const std::vector<offsets_container>& boundaries) const {
        std::vector<grains_container> res;
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

            res[num_grains - 1].push_back(off);
        }        
        return res;
    }
    grouped2_offsets_container sorted2_offsets() const {
        auto groupedbynum = group_boundaries_by_grains_num(boundaries_offsets());
        grouped2_offsets_container res;
        res.reserve(groupedbynum);
        for (auto& group : groupedbynum)
            res.emplace_back(std::move(group_boundaries_by_grains(group)));
        return res;
    }
    grouped2_grains_container grouped2_offsets_to_grains(const grouped2_offsets_container& bndoffsets) const {
        grouped2_grains_container res;
        res.reserve(bndoffsets.size());
        for (auto& samenum : bndoffsets) {
            res.emplace_back();
            res.back().reserve(samenum.size());
            for (auto& samebnd : samenum)
                res.back().push_back(get_grains(samebnd.front()));
        }
        return res;
    }

    std::vector<grains_container> supreme_pjoints_grains(const grouped2_grains_container& grconts) const {
        if (grconts.size() < 4)
            return {};

        std::vector<grains_container> res;
        for (auto samenum_it = grconts.begin() + 3; samenum_it < grconts.end(); ++samenum_it)
            for (auto& grains : *samenum_it)
                if (is_supreme_pjoint(grains, grconts))
                    res.push_back(grains);
        return res;
    }
    pos_type supreme_pjoint_pos(const offsets_container& offsets) const {
        std::size_t acc = 0;
        for (auto off : offsets)
            acc += off;
        return static_cast<Real>(acc) / offsets->size();
    }
    pos_type supreme_pjoint_pos(const grains_container& pjgrains, const grouped2_offsets_container& bndoffsets) const {
        auto& samenum = bndoffsets[pjgrains.size() - 1];
        auto it = grains_find_in_offsets(samenum, pjgrains);
        return supreme_pjoint_pos(*it);
    }
    std::vector<pos_type> supreme_pjoints_poses(
        const std::vector<grains_container>& pjgrains, const grouped2_offsets_container& bndoffsets) const {
        std::vector<pos_type> res;
        for (auto& pjgrs : pjgrains)
            res.push_back(supreme_pjoint_pos(pjgrs));
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
