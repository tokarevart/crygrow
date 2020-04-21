// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <optional>
#include "grain.h"
#include "neighborhood.h"
#include <unordered_set>


namespace cgr {

// clr means cellular
template <std::size_t Dim, typename Real = double>
class clr_grain {
public:
    using grain_type = cgr::grain<Dim, Real>;

    const grain_type* grain() const {
        return m_grain;
    }

    std::size_t range() const {
        return m_range;
    }
    void set_range(std::size_t range) {
        m_range = range;
        m_bbox_range = range * 2;
        m_shifts = nbh::make_shifts<Dim>(m_normfn, m_range, m_bbox_range);
    }

    const std::vector<std::size_t>& front() const {
        return m_front;
    }

    void front_swap_remove(std::size_t idx) {
        std::swap(m_front[idx], m_front.back());
        m_front.pop_back();
    }

    template <typename CrystedFn, typename NumGrainsFn>
    void advance_front(CrystedFn crysted, NumGrainsFn numgrs) {
        for (std::size_t i = 0; i < m_front.size();) {
            if (numgrs(m_front[i]) > 1)
                front_swap_remove(i);
            else
                ++i;
        }

        std::unordered_set<std::size_t> new_front;
        while (!m_front.empty()) {
            auto poses = apply_shifts(upos(front_pop()));
            for (auto& p : poses) {
                std::size_t o = offset(p);
                if (!crysted(o))
                    new_front.insert(o);
            }
        }
        m_front.assign(new_front.begin(), new_front.end());
    }

    template <typename InnerFn>
    void extract_front_from(std::vector<std::size_t> offs, InnerFn innfn, std::size_t thickness = 1) {
        m_front = std::move(offs);
        thin_front(innfn, thickness);
    }

    template <typename InnerFn>
    void thin_front(InnerFn innfn, std::size_t thickness = 1) {
        auto prev_front = std::move(m_front);
        auto shs = nbh::make_shifts<Dim>(norm_euclid<Dim>, thickness, thickness);
        for (std::size_t off : prev_front) {
            auto pos = static_cast<pos_t<Dim>>(upos(off));
            pos_t<Dim> del = pos_t<Dim>::filled_with(thickness);
            pos_t<Dim> corner0 = pos - del;
            pos_t<Dim> corner1 = pos + del;
            auto poses = inside(corner0) && inside(corner1) ?
                nbh::apply_shifts(pos, shs) :
                nbh::apply_shifts(pos, shs, [this](const pos_t<Dim>& p) -> bool { return inside(p); });
            bool near_boundary = false;
            for (auto& p : poses) {
                if (!innfn(offset(p), grain())) {
                    near_boundary = true;
                    break;
                }
            }
            if (near_boundary)
                m_front.push_back(off);
        }
    }

    clr_grain(const grain_type* grain, nbh::nbhood_kind kind, const upos_t<Dim>& dimlens, std::size_t nucleus_off)
        : m_grain{ grain }, m_dim_lens{ dimlens }, m_front{ nucleus_off } {
        switch (kind) {
        case nbh::nbhood_kind::von_neumann:
            m_normfn = norm_taxicab<Dim>; break;

        case nbh::nbhood_kind::moore:
            m_normfn = norm_chebyshev<Dim>; break;

        case nbh::nbhood_kind::euclid:
            m_normfn = norm_euclid<Dim>; break;

        case nbh::nbhood_kind::crystallographic:
            m_normfn = make_norm_cryst_fn<Dim, Real>(orientate_grow_dirs()); break;

        default:
            std::terminate();
        }
    }


private:
    const grain_type* m_grain;
    norm_fn<Dim> m_normfn;
    std::vector<pos_t<Dim>> m_shifts;
    std::size_t m_range = 0;
    std::size_t m_bbox_range = 0;
    upos_t<Dim> m_dim_lens;
    std::vector<std::size_t> m_front;

    upos_t<Dim> upos(std::size_t off) const {
        return cgr::upos(off, m_dim_lens);
    }
    std::size_t offset(const pos_t<Dim>& pos) const {
        return cgr::offset(pos, m_dim_lens);
    }

    bool inside(const pos_t<Dim>& pos) const {
        for (auto& e : pos.x)
            if (e < 0)
                return false;
        return cgr::inside(static_cast<upos_t<Dim>>(pos), m_dim_lens);
    }

    void front_push(std::size_t off) {
        m_front.push_back(off);
    }
    std::size_t front_pop() {
        std::size_t tmp = m_front.back();
        m_front.pop_back();
        return tmp;
    }

    std::vector<pos_t<Dim>> apply_shifts(const pos_t<Dim>& pos) const {
        pos_t<Dim> del = pos_t<Dim>::filled_with(m_bbox_range);
        pos_t<Dim> corner0 = pos - del;
        pos_t<Dim> corner1 = pos + del;
        if (inside(corner0) && inside(corner1))
            return nbh::apply_shifts(pos, m_shifts);
        else
            return nbh::apply_shifts(pos, m_shifts, [this](const pos_t<Dim>& pos) -> bool { return inside(pos); });
    }

    std::vector<grow_dir_t<Dim, Real>> orientate_grow_dirs() const {
        auto transposed_orien = m_grain->orientation().transposed();
        auto growdirs = m_grain->material()->grow_dirs();
        for (auto& gd : growdirs)
            gd = spt::dot(transposed_orien, gd);
        return growdirs;
    }
};

} // namespace cgr
