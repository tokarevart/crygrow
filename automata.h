// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <algorithm>
#include <map>
#include <cstddef>
#include <memory>
#include <optional>
#include <set>
#include <map>
#include <algorithm>
#include "neighborhood.h"
#include "vec.h"
#include "sptops.h"
#include "cgralgs.h"
#include "cell.h"
#include "clr-grain.h"


namespace cgr {

template <std::size_t Dim, typename Real = double>
class automata {
public:
    static constexpr std::size_t dim = Dim;
    using cell_type = cgr::cell<Dim, Real>;
    using grain_type = typename cell_type::grain_type;
    using grains_container = typename cell_type::grains_container;
    using material_type = typename grain_type::material_type;
    using orientation_type = typename grain_type::orientation_type;
    using clr_grain_type = clr_grain<Dim, Real>;
    using grow_dir_type = grow_dir_t<Dim, Real>;

    std::size_t num_crysted_cells() const {
        std::size_t res = 0;
        for (std::size_t i = 0; i < num_cells(); ++i)
            if (cell(i) && cell(i)->crysted)
                ++res;
        return res;
    }
    std::size_t num_cells() const {
        return m_cells.size();
    }
    const std::vector<cell_type>& cells() const {
        return m_cells;
    }

    const cell_type* cell(std::size_t offset) const {
        return m_cells[offset];
    }
    const cell_type* cell(const upos_t<Dim>& pos) const {
        return cell(offset(pos));
    }
    const cell_type* cell(const pos_t<Dim>& pos) const {
        return cell(offset(pos));
    }

    bool inside(const upos_t<Dim>& pos) const {
        return cgr::inside(pos, m_dim_lens);
    }
    bool inside(const pos_t<Dim>& pos) const {
        for (auto& e : pos.x)
            if (e < 0)
                return false;
        return cgr::inside(static_cast<upos_t<Dim>>(pos), m_dim_lens);
    }

    upos_t<Dim> upos(std::size_t offset) const {
        return cgr::upos(offset, m_dim_lens);
    }
    std::size_t offset(const upos_t<Dim>& pos) const {
        return offset(static_cast<pos_t<Dim>>(pos));
    }
    std::size_t offset(const pos_t<Dim>& pos) const {
        return cgr::offset(pos, m_dim_lens);
    }

    const upos_t<Dim>& dim_lens() const {
        return m_dim_lens;
    }
    std::size_t range() const {
        return m_range;
    }
    void set_range(std::size_t rng) {
        m_range = rng;
        for (auto& clrg : m_clrgrains)
            clrg.set_range(m_range);
    }

    bool stop_condition() const {
        for (std::size_t i = 0; i < num_cells(); ++i)
            if (!cell(i) || !cell(i)->crysted)
                return false;
        return true;
    }
    bool iterate() {
        if (stop_condition())
            return false;

        #pragma omp parallel for
        for (std::int64_t i = 0; i < m_clrgrains.size(); ++i)
            m_clrgrains[i].advance_front(
                [this](std::size_t off) -> bool { return m_cells[off] && m_cells[off]->crysted; },
                [this](std::size_t off) -> std::size_t { return !m_cells[off] ? 0 : m_cells[off]->grains.size(); });

        for (auto& clrg : m_clrgrains) {
            for (std::size_t off : clrg.front()) {
                if (!m_cells[off]) {
                    m_cells[off] = m_unicells[{ clrg.grain() }].get();
                } else {
                    std::set<const grain_type*> grs(m_cells[off]->grains.begin(), m_cells[off]->grains.end());
                    grs.insert(clrg.grain());
                    auto pcell = std::make_unique<cell_type>(nullptr, true);
                    pcell->grains.assign(m_cells[off]->grains.begin(), m_cells[off]->grains.end());
                    pcell->grains.push_back(clrg.grain());
                    auto [it, success] = m_unicells.insert({ grs, std::move(pcell) });
                    m_cells[off] = it->second.get();
                }
            }
        }

        for (auto& clrg : m_clrgrains)
            clrg.thin_front(
                [this](std::size_t off, const grain_type* gr) -> bool {
                    return m_cells[off] && m_cells[off]->grains.front() == gr && m_cells[off]->grains.size() == 1;
                });

        return true;
    }

    void spawn_grain(const grain_type* grain, const upos_t<Dim>& nucleus_pos, nbh::nbhood_kind kind) {
        spawn_grain(grain, offset(nucleus_pos), kind);
    }
    void spawn_grain(const grain_type* grain, std::size_t nucleus_off, nbh::nbhood_kind kind) {
        m_clrgrains.emplace_back(grain, kind, m_dim_lens, nucleus_off);
        m_clrgrains.back().set_range(m_range);
        auto [it, success] = m_unicells.insert({ std::set{ grain }, std::make_unique<cell_type>(grain, true) });
        m_cells[nucleus_off] = it->second.get();
    }

    void smooth(std::size_t rng) {
        auto shs = nbh::make_shifts<Dim>(norm_euclid<Dim>, rng);
        std::map<std::set<const grain_type*>, std::vector<std::size_t>> grconts;
        for (std::size_t i = 0; i < num_cells(); ++i) {
            auto pos = static_cast<pos_t<Dim>>(upos(i));
            std::vector<pos_t<Dim>> nbs = nbh::apply_shifts(pos, shs, 
                [this](const pos_t<Dim>& pos) -> bool { return inside(pos); });

            std::set<const grain_type*> grs(m_cells[i]->grains.begin(), m_cells[i]->grains.end());
            for (auto& nbpos : nbs) {
                auto& nbgrs = cell(nbpos)->grains;
                if (nbgrs.size() == 1)
                    grs.insert(nbgrs.front());
            }
            auto [it, success] = grconts.insert({ grs, std::vector<std::size_t>{} });
            it->second.push_back(i);
        }

        for (auto& p : grconts) {
            for (std::size_t i : p.second) {
                std::set<const grain_type*> grs(m_cells[i]->grains.begin(), m_cells[i]->grains.end());
                if (grs == p.first)
                    continue;

                auto pcell = std::make_unique<cell_type>(nullptr, true);
                pcell->grains.assign(p.first.begin(), p.first.end());
                auto [it, success] = m_unicells.insert({ p.first, std::move(pcell) });
                m_cells[i] = it->second.get();
            }
        }
    }

    void thin_boundary(std::size_t rng, std::size_t step = 1) {
        while (true) {
            if (rng + step >= range())
                set_range(rng);
            else
                set_range(range() - step);
            std::cout << "started range=" << range() << std::endl;

            #pragma omp parallel for
            for (std::int64_t i = 0; i < m_clrgrains.size(); ++i) {
                std::vector<std::size_t> inner_cells;
                for (std::size_t j = 0; j < num_cells(); ++j)
                    if (m_cells[j]->grains.front() == m_clrgrains[i].grain() &&
                        m_cells[j]->grains.size() == 1)
                        inner_cells.push_back(j);

                m_clrgrains[i].extract_front_from(inner_cells,
                    [this](std::size_t off, const grain_type* gr) -> bool {
                        return m_cells[off]->grains.front() == gr && m_cells[off]->grains.size() == 1;
                    });
            }

            for (auto& cl : m_cells)
                if (cl->grains.size() > 1)
                    cl = nullptr;

            while (!stop_condition()) {
                iterate();
                bool all_empty_fronts = true;
                for (auto& clrgr : m_clrgrains) {
                    if (!clrgr.front().empty()) {
                        all_empty_fronts = false;
                        break;
                    }
                }
                if (all_empty_fronts)
                    while (!extrapolate_nullcells());
            }

            std::size_t num_single_grained_cells = 0;
            for (const cell_type* cl : m_cells)
                if (cl->grains.size() == 1)
                    ++num_single_grained_cells;
            if (range() == rng || num_single_grained_cells == num_cells())
                break;
        }
        while (!extrapolate_cells_with_numgrains_gt2());
    }

    Real diam(std::vector<pos_t<Dim>> poses) const {
        std::int64_t max_diff2 = 0.0;
        for (pos_t<Dim> p0 : poses) {
            std::int64_t maxd2 = 0;
            for (pos_t<Dim> p1 : poses) {
                std::int64_t d2 = (p1 - p0).magnitude2();
                maxd2 = std::max(maxd2, d2);
            }
            max_diff2 = std::max(max_diff2, maxd2);
        }

        return std::sqrt(max_diff2);
    }

    Real cells_diam(cell_type* cell) const {
        std::vector<pos_t<Dim>> poses;
        for (std::size_t i = 0; i < m_cells.size(); ++i)
            if (cell == m_cells[i])
                poses.push_back(static_cast<pos_t<Dim>>(upos(i)));

        return diam(poses);
    }

    std::vector<Real> diams_inter4() const {
        std::vector<Real> res;
        for (auto& [grs, upcell] : m_unicells) {
            if (grs.size() != 4)
                continue;

            res.push_back(cells_diam(upcell.get()));
        }

        return res;
    }

    automata(std::size_t dimlen)
        : automata(upos_t<Dim>::filled_with(dimlen)) {}
    automata(const upos_t<Dim>& dimlens) : m_dim_lens{ dimlens } {
        std::size_t new_num_cells = std::accumulate(
            m_dim_lens.x.begin(), m_dim_lens.x.end(),
            static_cast<std::size_t>(1), std::multiplies<std::size_t>());
        m_cells.assign(new_num_cells, nullptr);
    }


private:
    std::size_t m_range = 0;
    upos_t<Dim> m_dim_lens;

    std::vector<cell_type*> m_cells;
    std::vector<clr_grain_type> m_clrgrains;
    std::map<std::set<const grain_type*>, std::unique_ptr<cell_type>> m_unicells;

    template <typename Pred>
    bool extrapolate_cells(Pred pred) {
        bool success = true;

        auto shs = nbh::make_shifts<Dim>(norm_euclid<Dim>, 1);
        for (std::size_t i = 0; i < num_cells(); ++i) {
            if (!pred(i)) continue;

            auto pos = static_cast<pos_t<Dim>>(upos(i));
            std::vector<pos_t<Dim>> nbs = nbh::apply_shifts(pos, shs,
                [this](const pos_t<Dim>& pos) -> bool { return inside(pos); });

            std::vector<std::pair<const cell_type*, std::size_t>> prs;
            auto prs_find_cell = [&prs](const cell_type* cl) -> std::size_t {
                std::size_t i = 0;
                for (; i < prs.size(); ++i)
                    if (prs[i].first == cl)
                        break;
                return i;
            };

            for (auto& nb : nbs) {
                const cell_type* nbcl = cell(nb);
                if (!nbcl || pred(offset(nb))) continue;

                std::size_t find_res = prs_find_cell(nbcl);
                if (find_res == prs.size())
                    prs.push_back({ nbcl, 1 });
                else
                    ++prs[find_res].second;
            }
            std::sort(prs.begin(), prs.end(),
                [](const auto& pr0, const auto& pr1) -> bool {
                    return pr0.second < pr1.second;
                });

            if (prs.empty()) {
                success = false;
                continue;
            }
            m_cells[i] = const_cast<cell_type*>(prs.back().first);
        }

        return success;
    }

    bool extrapolate_cells_with_numgrains_gt2() {
        return extrapolate_cells(
            [this](std::size_t off) -> bool {
                return m_cells[off] && m_cells[off]->grains.size() > 2;
            });
    }

    bool extrapolate_nullcells() {
        return extrapolate_cells(
            [this](std::size_t off) -> bool {
                return !m_cells[off];
            });
    }
};

} // namespace cgr
