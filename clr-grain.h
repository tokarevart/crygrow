// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <optional>
#include "grain.h"
#include "neighborhood.h"


namespace cgr {

// clr means cellular
template <std::size_t Dim, typename Real = double>
class clr_grain {
public:
    using grain_type = cgr::grain<Dim, Real>;
    using crysted_fn = std::function<bool(std::size_t)>;

    const grain_type* grain() const {
        return m_grain;
    }

    std::size_t range() const {
        return m_range;
    }
    void set_range(std::size_t range) {
        m_shifts = make_shifts<Dim>(m_normfn, range);
    }

    const std::vector<std::size_t>& front() const {
        return m_front;
    }

    clr_grain(const grain_type* grain, nbh::nbhood_kind kind, const upos_t<Dim>& dimlens, std::size_t nucleus_off, crysted_fn crfn)
        : m_grain{ grain }, m_orientation{ orien }, m_dim_lens{ dimlens }, m_front{ nucleus_off }, m_crfn{ crfn } {
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
    upos_t<Dim> m_dim_lens;
    crysted_fn m_crfn;
    std::vector<std::size_t> m_front;

    bool crysted(std::size_t off) const {
        return m_crfn(off);
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
        return nbh::apply_shifts(pos, m_shifts, inside);
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
