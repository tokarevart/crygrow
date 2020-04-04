// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "grain.h"
#include "neighborhood.h"

namespace cgr {

template <std::size_t Dim, typename Real = double>
class clr_grain {
public:
    using grain_type = cgr::grain<Dim, Real>;
    using orientation_type = orientation_t<Dim, Real>;

    const grain_type* material() const {
        return m_material;
    }

    clr_grain(const grain_type* grain, nbhood_kind kind)
        : m_grain{ grain }, m_orientation{ orien } {
        switch (kind) {
        case nbhood_kind::von_neumann:
            m_normfn = norm_taxicab; break;

        case nbhood_kind::moore:
            m_normfn = norm_chebyshev; break;

        case nbhood_kind::euclid:
            m_normfn = norm_euclid; break;

        case nbhood_kind::crystallographic:
            m_normfn = make_norm_cryst_fn(orientate_grow_dirs()); break;

        default:
            std::terminate();
        }
    }


private:
    const grain_type* m_grain;
    norm_fn<Dim> m_normfn;

    std::vector<spt::vec<Dim, Real>> orientate_grow_dirs() const {
        auto transposed_orien = m_grain->orientation().transposed();
        auto growdirs = m_grain->material()->grow_dirs();
        for (auto& gd : growdirs)
            gd = spt::dot(transposed_orien, gd);
        return growdirs;
    }
};

} // namespace cgr
