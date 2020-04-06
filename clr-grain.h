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
    using offdel_pair = std::pair<std::size_t, Real>; // offset and delta
    using offset_fn = std::function<std::size_t(const spt::veci<Dim>&)>;
    using crysted_fn = std::function<bool(std::size_t)>;

    const grain_type* grain() const {
        return m_grain;
    }

    void range() const {
        return m_range;
    }
    void set_range(std::size_t range) {
        m_shifts = make_shifts<Dim>(m_normfn, range);
    }

    void set_insidefn(inside_fn infn) {
        m_infn = infn;
    }
    void set_offsetfn(offset_fn offfn) {
        m_offfn = offfn;
    }

    offdel_pair next_offdel() {
        offdel_pair temp = m_offdels.back();
        m_offdels.pop_back();
        return temp;
    }
    std::optional<offdel_pair> try_next_offdel() {
        if (empty_offdels())
            return std::nullopt;
        else
            return next_offdel();
    }
    bool empty_offdels() const {
        return m_offdels.empty();
    }

    clr_grain(const grain_type* grain, nbhood_kind kind)
        : m_grain{ grain }, m_orientation{ orien } {
        switch (kind) {
        case nbhood_kind::von_neumann:
            m_normfn = norm_taxicab<Dim>; break;

        case nbhood_kind::moore:
            m_normfn = norm_chebyshev<Dim>; break;

        case nbhood_kind::euclid:
            m_normfn = norm_euclid<Dim>; break;

        case nbhood_kind::crystallographic:
            m_normfn = make_norm_cryst_fn<Dim>(orientate_grow_dirs()); break;

        default:
            std::terminate();
        }
    }


private:
    const grain_type* m_grain;
    norm_fn<Dim> m_normfn;
    nbh::poses_t<Dim> m_shifts;
    std::size_t m_range = 0;
    nbh::inside_fn m_infn;
    offset_fn m_offfn;
    crysted_fn m_crfn;
    std::vector<offdel_pair> m_offdels;

    bool crysted(std::size_t off) const {
        return m_crfn(off);
    }

    void add_offdel(offdel_pair od) {
        m_offdels.push_back(od);
    }

    nbh::poses_t<Dim> apply_shifts()

    std::vector<spt::vec<Dim, Real>> orientate_grow_dirs() const {
        auto transposed_orien = m_grain->orientation().transposed();
        auto growdirs = m_grain->material()->grow_dirs();
        for (auto& gd : growdirs)
            gd = spt::dot(transposed_orien, gd);
        return growdirs;
    }
};

} // namespace cgr
