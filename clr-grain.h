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
    using orientation_type = orientation_t<Dim, Real>;
    using offgrain_pair = std::pair<std::size_t, const grain_type*>;

    const grain_type* grain() const {
        return m_grain;
    }

    void range() const {
        return m_range;
    }
    void set_range(std::size_t range) {
        m_shifts = make_nbhood_pos_shifts<Dim>(m_normfn, range);
    }

    offgrain_pair next_offgrain() {
        offgrain_pair temp = m_offgrains.back();
        m_offgrains.pop_back();
        return temp;
    }
    std::optional<offgrain_pair> try_next_offgrain() {
        if (empty_offgrains())
            return std::nullopt;
        else
            return next_offgrain();
    }
    bool empty_offgrains() const {
        return m_offgrains.empty();
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
    nbhood_pos<Dim> m_shifts;
    std::size_t m_range = 0;
    std::vector<offgrain_pair> m_offgrains;

    std::vector<spt::vec<Dim, Real>> orientate_grow_dirs() const {
        auto transposed_orien = m_grain->orientation().transposed();
        auto growdirs = m_grain->material()->grow_dirs();
        for (auto& gd : growdirs)
            gd = spt::dot(transposed_orien, gd);
        return growdirs;
    }
};

} // namespace cgr
