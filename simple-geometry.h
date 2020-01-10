// Copyright © 2019-2020 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <optional>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include "simple-automata.h"
#include "grgeo.h"


namespace cgr {

class simple_geo_tool {
public:
    static constexpr std::size_t dim = 3;
    using vecu = spt::vec3u;
    using real_type = grgeo::real_type;
    using tag_type = grgeo::tag_type;
    using utag_type = grgeo::utag_type;
    using pos_type = grgeo::pos_type;

    using cell_type = grgeo::cell_type;
    using cells_container = grgeo::cells_container;
    using offsets_container = std::vector<std::size_t>;
    using grain_type = cell_type::grain_type;
    using grains_container = cell_type::grains_container;
    using vertices_container = offsets_container;
    using edges_container = std::vector<offsets_container>;
    using faces_container = std::vector<offsets_container>;
    // grains num then boundaries grouped
    template <typename T>
    using vector2gd = std::vector<std::vector<T>>;

    using automata_type = cgr::simple_automata<dim, nbhood_kind::euclid, real_type>;
    using gr_geometry = grgeo::gr_geometry;

    void add_empty_gr_volumes(const std::vector<grains_container>& grconts) {
        std::unordered_set<grain_type*> uniquegrs;
        for (auto& grcont : grconts)
            for (auto& gr : grcont)
                if (gr->material())
                    uniquegrs.insert(gr);

        for (auto it = uniquegrs.begin(); it != uniquegrs.end(); ++it)
            m_gr_geo.add_gr_volume(*it);
    }
    void add_empty_gr_surfaces(const std::vector<grains_container>& grconts) {
        for (auto& grcont : grconts)
            m_gr_geo.add_gr_surface(&grcont);
    }
    void add_empty_gr_lines(const std::vector<grains_container>& grconts) {
        for (auto& grcont : grconts)
            m_gr_geo.add_gr_line(&grcont);
    }
    void add_gr_points(const std::vector<grains_container>& grconts,
                       const std::vector<offsets_container>& offconts) {
        for (std::size_t i = 0; i < grconts.size(); ++i)
            m_gr_geo.add_gr_point(&grconts[i], central_pos(offconts[i]));
    }
    
    bool grains_contains_unsorted(const grains_container& first, const grain_type* second) const {
        return std::find(first.begin(), first.end(), second) != first.end();
    }
    bool grains_includes_unsorted(const grains_container& first, const grains_container& second) const {
        if (first.size() < second.size())
            return false;
        for (auto s : second)
            if (!grains_contains_unsorted(first, s))
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
        std::vector<offsets_container>& offconts, const grains_container& grains) const {
        auto it = offconts.begin();
        for (; it < offconts.end(); ++it)
            if (grains_equal_unsorted(get_grains(it->front()), grains))
                break;
        return it;
    }
    
    template <std::size_t Dir>
    std::size_t shift(std::size_t origin, std::size_t dist) const {
        if constexpr (Dir == 0) 
            return origin + dist;

        if constexpr (Dir == 1) 
            return origin + get_dim_lens()[0] * dist;

        if constexpr (Dir == 2) 
            return origin + get_dim_lens()[1] * get_dim_lens()[0] * dist;
    }
    // exclusive dist_end
    template <std::size_t Dir>
    offsets_container shift(std::size_t origin, std::size_t dist_begin, std::size_t dist_end) const {
        offsets_container res;
        res.reserve(dist_end - dist_begin);
        for (std::size_t i = dist_begin; i < dist_end; ++i)
            res.push_back(shift<Dir>(origin, i));
        return res;
    }
    template <std::size_t Dir>
    offsets_container shift(offsets_container origins, std::size_t dist_begin, std::size_t dist_end) const {
        offsets_container res;
        res.reserve((dist_end - dist_begin) * origins.size());
        for (auto origin : origins)
            for (std::size_t i = dist_begin; i < dist_end; ++i)
                res.push_back(shift<Dir>(origin, i));
        return res;
    }

    template <std::size_t Dir>
    std::size_t vshift(std::size_t origin) const {
        return shift<Dir>(origin, get_dim_lens()[Dir] - 1);
    }
    template <std::size_t Dir>
    offsets_container eshift(std::size_t vertex) const {
        return shift<Dir>(vertex, 0, get_dim_lens()[Dir]);
    }
    template <std::size_t Dir>
    offsets_container fshift(const offsets_container& edge) const {
        return shift<Dir>(edge, 0, get_dim_lens()[Dir]);
    }

    vertices_container box_vertices() const {
        vertices_container v(8);
        v[0] = 0;
        v[1] = vshift<0>(v[0]);
        v[2] = vshift<1>(v[0]);
        v[3] = vshift<0>(v[2]);
        v[4] = vshift<2>(v[0]);
        v[5] = vshift<0>(v[4]);
        v[6] = vshift<1>(v[4]);
        v[7] = vshift<0>(v[6]);
        return v;
    }
    edges_container box_edges(const vertices_container& v) const {
        edges_container e(12);
        e[0] = std::move(eshift<0>(v[0]));
        e[1] = std::move(eshift<1>(v[0]));
        e[2] = std::move(eshift<0>(v[2]));
        e[3] = std::move(eshift<1>(v[1]));
        e[4] = std::move(eshift<0>(v[4]));
        e[5] = std::move(eshift<1>(v[4]));
        e[6] = std::move(eshift<0>(v[6]));
        e[7] = std::move(eshift<1>(v[5]));
        e[8] = std::move(eshift<2>(v[0]));
        e[9] = std::move(eshift<2>(v[1]));
        e[10] = std::move(eshift<2>(v[2]));
        e[11] = std::move(eshift<2>(v[3]));
        return e;
    }
    faces_container box_faces(const edges_container& e) const {
        faces_container f(6);
        f[0] = std::move(fshift<2>(e[0]));
        f[1] = std::move(fshift<2>(e[1]));
        f[2] = std::move(fshift<2>(e[2]));
        f[3] = std::move(fshift<2>(e[3]));
        f[4] = std::move(fshift<1>(e[4]));
        f[5] = std::move(fshift<1>(e[0]));
        return f;
    }

    offsets_container flatten(const std::vector<offsets_container>& conts) const {
        offsets_container res(std::accumulate<std::size_t>(
            conts.begin(), conts.end(), 0,
            [](const offsets_container& e) { return e.size(); }));
        for (auto& cont : conts)
            for (auto off : cont)
                res.push_back(off);
        return res;
    }
    offsets_container flatten_unique(const std::vector<offsets_container>& conts) const {
        std::unordered_set<std::size_t> us;
        for (auto& cont : conts)
            for (auto off : cont)
                us.insert(off);
        return offsets_container(us.begin(), us.end());
    }
    
    offsets_container boundaries_offsets() const {
        offsets_container res;
        for (std::size_t i = 0; i < num_cells(); ++i) {
            if (get_grains(i).size() > 1)
                res.push_back(i);
            // dirty hack
            if (get_grains(i).size() > 4)
                std::cout << "oops!";
        }
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
    std::vector<offsets_container> group_boundaries_by_grains(const offsets_container& bryoffsets) const {
        std::vector<offsets_container> res;
        for (auto off : bryoffsets) {
            auto found = grains_find_in_offsets(res, get_grains(off));
            if (found == res.end()) {
                res.push_back(offsets_container{ off });
            } else {
                found->push_back(off);
            }
        }
        return res;
    }
    std::vector<offsets_container> group_boundaries_by_grains_num(const offsets_container& bryoffsets) const {
        std::vector<offsets_container> res;
        for (auto off : bryoffsets) {
            std::size_t num_grains = get_grains(off).size();
            if (num_grains > res.size())
                while (res.size() < num_grains - 1)
                    res.emplace_back();

            res[num_grains - 2].push_back(off);
        }        
        return res;
    }
    vector2gd<offsets_container> grouped2_offsets() const {
        auto groupedbynum = group_boundaries_by_grains_num(boundaries_offsets());
        vector2gd<offsets_container> res;
        res.reserve(groupedbynum.size());
        for (auto& group : groupedbynum)
            res.push_back(std::move(group_boundaries_by_grains(group)));
        return res;
    }
    vector2gd<grains_container> grouped2_offsets_to_grains(const vector2gd<offsets_container>& offsconts) const {
        vector2gd<grains_container> res;
        res.reserve(offsconts.size());
        for (auto& samenum : offsconts)
            res.push_back(std::move(boundaries_grains(samenum)));
        return res;
    }

    pos_type central_pos(const offsets_container& offsets) const {
        vecu acc;
        for (auto off : offsets)
            acc += cgr::upos(off, get_dim_lens());
        auto accreal = static_cast<pos_type>(acc);
        for (auto& e : accreal.x)
            e /= offsets.size();
        return accreal;
    }


    void connect_gr_geometry() {
        for (auto& vol : m_gr_geo.gr_volumes)
            for (auto& face : m_gr_geo.gr_surfaces)
                if (grains_contains_unsorted(*face->pgrains, vol->pgrain))
                    vol->faces.emplace_back(std::make_unique<oriented<gr_face>>(face.get()));

        for (auto& face : m_gr_geo.gr_surfaces)
            for (auto& edge : m_gr_geo.gr_lines)
                if (grains_includes_unsorted(*edge->pgrains, *face->pgrains))
                    face->edges.emplace_back(std::make_unique<oriented<gr_edge>>(edge.get()));

        for (auto& edge : m_gr_geo.gr_lines)
            for (auto& vert : m_gr_geo.gr_points)
                if (grains_includes_unsorted(*vert->pgrains, *edge->pgrains))
                    edge->verts.push_back(vert.get());
    }

    void init_gr_geometry(const std::vector<offsets_container>& pjointsoffs) {
        add_empty_gr_volumes(m_g2grs[0]);
        add_empty_gr_surfaces(m_g2grs[0]);
        add_empty_gr_lines(m_g2grs[1]);
        add_gr_points(m_g2grs[2], pjointsoffs);
        connect_gr_geometry();
    }

    void make_geometry() {
        auto box_vs = box_vertices();
        auto box_es = box_edges(box_vs);
        auto box_fs = box_faces(box_es);

        for (std::size_t i = 0; i < 6; ++i) {
            m_boxbry_grains[i] = std::make_unique<grain_type>(nullptr);
            for (auto off : box_fs[i]) {
                add_grain(off, m_boxbry_grains[i].get());
            }
        }

        auto g2offs = grouped2_offsets();
        m_g2grs = grouped2_offsets_to_grains(g2offs);
        init_gr_geometry(g2offs[2]);
        g2offs.clear();
        g2offs.shrink_to_fit();

        m_gr_geo.geometry.orient_lines();
        //
    }
    
    void write_geo(std::ostream& os) const {
        m_gr_geo.geometry.write(os);
    }


    std::optional<std::string> is_global_max_order_overflow() const {
        for (std::size_t i = 0; i < num_cells(); ++i)
            if (get_grains(i).size() > 4)
                return "Max boundary order overflow\n";
        return std::nullopt;
    }
    std::optional<std::string> is_box_vertices_max_order_overflow() const {
        //
    }
    std::optional<std::string> is_box_edges_max_order_overflow() const {
        //
    }
    std::optional<std::string> is_box_faces_max_order_overflow() const {
        //
    }
    std::optional<std::string> is_exception() const {
        //
    }

    simple_geo_tool(const automata_type* automata)
        : m_automata(automata) {}


private:
    const automata_type* m_automata;

    std::unordered_map<std::size_t, std::unique_ptr<grains_container>> m_boxbry_grconts;
    std::array<std::unique_ptr<grain_type>, 6> m_boxbry_grains;
    vector2gd<grains_container> m_g2grs;
    gr_geometry m_gr_geo;

    std::size_t num_cells() const {
        return m_automata->num_cells();
    }
    cell_type* get_cell(std::size_t offset) const {
        return m_automata->get_cell(offset);
    }
    void add_grain(std::size_t offset, const grain_type* pgrain) {
        m_boxbry_grconts.insert({ offset, std::make_unique<grains_container>(get_cell(offset)->grains) });
        m_boxbry_grconts[offset]->push_back(const_cast<grain_type*>(pgrain));
    }
    const grains_container& get_grains(std::size_t offset) const {
        vecu dlens = get_dim_lens();
        spt::vec3u upos = cgr::upos(offset, dlens);

        for (std::size_t i = 0; i < 3; ++i)
            if (upos[i] == 0 ||
                upos[i] == dlens[i] - 1) {
                auto searchres = m_boxbry_grconts.find(offset);
                return const_cast<const grains_container&>(*searchres->second);
            }

        return get_cell(offset)->grains;
    }
    const vecu& get_dim_lens() const {
        return m_automata->get_dim_lens();
    }
};

} // namespace cgr
