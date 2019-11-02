#pragma once
#include <algorithm>
#include <set>
#include "automata-base.h"
#include "simple-cell.h"


namespace cgr {

template <std::size_t Dim, nbhood_kind NbhoodKind = nbhood_kind::euclid, typename Real = double>
class simple_automata
    : public automata_base<Dim, NbhoodKind, simple_cell<Dim, Real>, cell_mut_group::mutable_only> {
public:
    static constexpr Real epsilon = std::numeric_limits<Real>::epsilon();
    using base = automata_base<Dim, NbhoodKind, simple_cell<Dim, Real>, cell_mut_group::mutable_only>;
    using veci = typename base::veci;
    using vecu = typename base::vecu;
    using cell_type = typename base::cell_type;
    using nbhood_type = typename base::nbhood_type;
    using nbhood_offset_type = typename base::nbhood_offset_type;
    using nbhood_offsets_container = typename base::nbhood_offsets_container;
    using cells_delta_container = std::vector<Real>;
    using crystallite_type = typename cell_type::crystallite_type;
    using material_type = typename crystallite_type::material_type;
    using orientation_type = typename crystallite_type::orientation_type;
    using grow_dir = typename material_type::grow_dir;

    bool stop_condition() const override {
        for (std::size_t i = 0; i < base::num_cells(); ++i)
            if (base::get_cell(i)->crystallinity < 1.0 - epsilon * (1.0 + base::get_cell(i)->crystallinity))
                return false;

        return true;
    }
    bool iterate() override {
        if (stop_condition())
            return false;

        if (!is_nbhood_offsets_initialized)
            initialize_nbhood_offsets();

        for (auto& delta : m_cells_delta)
            delta = 0.0;

        auto defnbhoffset = base::make_nbhood_offset();
        Real accdefdpmagn = 0.0;
        for (auto nboff : defnbhoffset) 
            accdefdpmagn += static_cast<Real>((
                veci::filled_with(base::default_range()) 
                - cgr::upos(nboff, vecu::filled_with(2 * base::default_range() + 1))
                ).magnitude());
        
        #pragma omp parallel 
        {
            #pragma omp for
            for (std::int64_t i = 0; i < static_cast<std::int64_t>(base::num_cells()); ++i) {
                auto curpos = base::pos(i);
                auto pcell = base::get_cell(i);
                if (std::abs(pcell->crystallinity - 1.0) <= epsilon * (1.0 + pcell->crystallinity))
                    continue;
                
                Real curdelta = 0.0;
                std::set<std::unique_ptr<grow_dir>> gds;
                grow_dir accgd;
                grow_dir accdp;
                Real accdpmagn = 0.0;
                Real accabsdot = 0.0;
                std::vector<grow_dir*> bufgds;
                for (auto nboff : base::get_nbhood_offset(i)) {
                    auto pnb = base::get_cell(nboff);
                    
                    if (pnb->crystallinity < 1.0 - epsilon * (1.0 + pcell->crystallinity) ||
                        pnb->crystallites.size() != 1) {
                        continue;
                    }

                    if (std::find(pcell->crystallites.begin(),
                                  pcell->crystallites.end(),
                                  pnb->crystallites.front()) == pcell->crystallites.end())
                        pcell->crystallites.push_back(pnb->crystallites.front());
                    
                    grow_dir deltapos = curpos - base::pos(nboff);
                    accdp += deltapos;
                    accdpmagn += deltapos.magnitude();

                    if (pnb->crystallites.front()->material()->matproperty() == material_property::anisotropic) {
                        auto tranorient = pnb->crystallites.front()->orientation().transposed();

                        for (auto& gd : pnb->crystallites.front()->material()->grow_dirs()) {
                            bufgds.emplace_back(new grow_dir(spt::dot(tranorient, gd)));
                            accabsdot += std::abs(spt::dot(deltapos, *bufgds.back()));
                        }
                        for (auto& gd : bufgds)
                            gds.insert(std::move(std::unique_ptr<grow_dir>(gd)));

                        for (auto& gd : bufgds)
                            accabsdot += std::abs(spt::dot(deltapos, *gd));

                        bufgds.clear();
                    } else {
                        gds.insert(std::make_unique<grow_dir>((deltapos).normalize()));
                        accabsdot += std::abs(spt::dot(deltapos, deltapos));
                    }
                }
                for (auto& gd : gds) {
                    curdelta += std::abs(spt::dot(accdp, *gd)) * accdp.magnitude() / accdpmagn;
                    Real f2_coef = pcell->crystallites.size() > 1 ? 3 : 1;
                    Real factor2 = std::clamp((f2_coef * accdpmagn - accdefdpmagn / 2) / (accdefdpmagn / 2), 0.0, 1.0);
                    curdelta += accabsdot * factor2;
                }
                
                if (curdelta > epsilon && !gds.empty()) 
                    curdelta /= gds.size() * 8;

                if (curdelta > epsilon)
                    m_cells_delta[i] += curdelta / base::default_nbhood_size();
            }

            #pragma omp barrier
            #pragma omp for
            for (std::int64_t i = 0; i < static_cast<std::int64_t>(base::num_cells()); ++i) {
                auto pcell = base::get_cell(i);
                pcell->crystallinity += m_cells_delta[i];
                m_cells_delta[i] = 0.0;

                if (pcell->crystallinity > 1.0 + epsilon * (1.0 + pcell->crystallinity)) {
                    pcell->crystallinity = 1.0;
                    update_nbhood_offsets_after_crystallization(i);
                }
            }
        }

        return true;
    }

    simple_automata(std::size_t dim_len,
                      std::size_t default_range = 1)
        : simple_automata(veci::zeros(), vecu::filled_with(dim_len), default_range) {}

    simple_automata(const vecu& dim_lens,
                      std::size_t default_range = 1)
        : simple_automata(veci::zeros(), dim_lens, default_range) {}

    simple_automata(const veci& corner0, const veci& corner1, 
                      std::size_t default_range = 1)
        : base(corner0, corner1, default_range) {
        m_cells_delta.assign(base::num_cells(), 0.0);
    }


private:
    cells_delta_container m_cells_delta;
    // todo: store common cells, such as empty cell and crystallized cell for each crystallite

    bool is_nbhood_offsets_initialized = false;
    void initialize_nbhood_offsets() {
        #pragma omp parallel for
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(base::num_cells()); ++i)
            if (std::abs(base::get_cell(i)->crystallinity - 1.0)
                <= epsilon * (1.0 + base::get_cell(i)->crystallinity)) {
                for (auto nboff : base::make_nbhood_offset(i)) {
                    if (base::get_cell(nboff)->crystallinity <= epsilon &&
                        base::get_nbhood_offset(nboff).empty())
                        base::set_nbhood_offset(nboff);
                }
            }

        is_nbhood_offsets_initialized = true;
    }
    void update_nbhood_offsets_after_crystallization(std::size_t offset) {
        for (auto nboff : base::get_nbhood_offset(offset)) {
            if (base::get_cell(nboff)->crystallinity <= epsilon &&
                base::get_nbhood_offset(nboff).empty())
                base::set_nbhood_offset(nboff);
        }

        if (!base::get_nbhood_offset(offset).empty())
            base::clear_nbhood_offset(offset);
    }
};

} // namespace cgr
