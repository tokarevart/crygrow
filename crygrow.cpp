// Copyright © 2019-2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#include <iostream>
#include <fstream>
#include <random>
#include "sptalgs.h"
#include "automata.h"
#include "geometry.h"
#include "progress-bar.h"

#define DIM3

#ifdef DIM3
constexpr std::size_t dim = 3;
#else
constexpr std::size_t dim = 2;
#endif
std::size_t seed = 0;
using automata_t = cgr::automata<dim>;
using cell_t = cgr::cell<dim>;
using grain_t = cgr::grain<dim>;
using material_t = cgr::material<dim>;


std::vector<cgr::upos_t<dim>> make_central_pos(std::size_t size) {
    std::vector<cgr::upos_t<dim>> res;
    res.push_back(cgr::upos_t<dim>::filled_with(size / 2));
    return res;
}

std::uint64_t min_distance2(const cgr::upos_t<dim>& pos, const std::vector<cgr::upos_t<dim>>& others) {
    std::uint64_t res = std::numeric_limits<std::uint64_t>::max();
    for (auto& other : others) {
        std::uint64_t dist = (other - pos).magnitude2();
        if (dist < res)
            res = dist;
    }
    return res;
}

std::vector<cgr::upos_t<dim>> make_random_poses(std::size_t size, std::size_t num, std::uint64_t min_dist2 = 0) {
    std::vector<cgr::upos_t<dim>> res;
    std::mt19937_64 gen(seed);
    std::uniform_int_distribution<std::size_t> dis(0, size - 1);
    for (std::size_t i = 0; i < num;) {
        cgr::upos_t<dim> curpos;
        for (auto& e : curpos.x)
            e = dis(gen);

        if (min_distance2(curpos, res) >= min_dist2) {
            bool flag = false;
            for (auto e : curpos.x)
                if (!(static_cast<std::size_t>(e * e) >= min_dist2 &&
                      static_cast<std::size_t>(e) < size - static_cast<std::size_t>(std::sqrt(min_dist2)))) {
                    flag = true;
                    break;
                }
            if (flag)
                continue;

            res.push_back(curpos);
            ++i;
        }
    }
    return res;
}

#include "range-iter.h"
int main_test() {
    for (itr::range_iter<int> it(0, 7); !it.has_ended(); ++it)
        std::cout << *it << ' ';
    std::cout << std::endl;
    for (itr::range_iter<int> it(0, 7, itr::dir::reverse); !it.has_ended(); ++it)
        std::cout << *it << ' ';

    return 0;
}

void write_image_pixels(std::ostream& os, const automata_t& atmt, bool blackwhite = false) {
    for (std::size_t i = 0; i < atmt.dim_lens()[0] * atmt.dim_lens()[1]; ++i) {
        auto curpos = atmt.upos(i);
        #ifdef DIM3
        //curpos[2] = static_cast<std::int64_t>(size) / 2;
        //curpos[2] = size - 1;
        curpos[2] = 0;
        #endif
        auto& cell = *atmt.cell(curpos);

        std::array<std::uint8_t, 3> color;
        if (!blackwhite) {
            if (!cell.crysted &&
                !cell.grains.empty()) {
                color = { 0, 255, 0 };
            } else if (!cell.crysted ||
                cell.grains.empty()) {
                color = { 255, 255, 255 };
            } else if (cell.grains.size() == 1) {
                color = { 0, 0, 0 };
            } else if (cell.grains.size() == 2) {
                color = { 0, 0, 255 };
            } else {
                color = { 255, 0, 0 };
            }
        } else {
            if (cell.grains.size() == 1 &&
                cell.crysted) {
                color = { 0, 0, 0 };
            } else {
                color = { 255, 255, 255 };
            }
        }


        os  << curpos[0] << ' '
            << curpos[1] << ' '
            << static_cast<int>(color[0]) << ' '
            << static_cast<int>(color[1]) << ' '
            << static_cast<int>(color[2]) << std::endl;
    }
}

int inner_main() {
    std::size_t size = 300;
    std::size_t range = 8;
    automata_t automata(size);
    automata.set_range(range);

    //auto init_poses = make_central_pos(size);
    auto init_poses = make_random_poses(size, 10, std::pow(range * 1, 2));

    //material_t mater;
    material_t mater({ 
        #ifdef DIM3
        spt::vecd<dim>({ 1.0, 0.0, 0.0 }).normalize(),
        spt::vecd<dim>({ 0.0, 1.0, 0.0 }).normalize(),
        spt::vecd<dim>({ 0.0, 0.0, 1.0 }).normalize() });
        #else
        spt::vecd<dim>({ 4.0, 1.0 }).normalize(), 
        spt::vecd<dim>({ -1.0, 4.0 }).normalize() });
        #endif
    //std::vector<grain_t> grains(init_poses.size(), grain_t(&mater));
    std::vector<grain_t> grains;
    grains.reserve(init_poses.size());
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> dis(-1.0, std::nextafter(1.0, 2.0));
    for (std::size_t i = 0; i < init_poses.size(); ++i) {
        #ifdef DIM3
        const double pi = 3.14159265359;
        auto rot = spt::rotation(spt::vecd<dim>({ dis(gen), dis(gen), dis(gen) }).normalize(), std::abs(dis(gen)) * pi);
        grains.emplace_back(&mater, rot);
        #else
        auto first = spt::vecd<dim>({ dis(gen), dis(gen) }).normalize();
        spt::vecd<dim> second({ -first[1], first[0] });
        grains.emplace_back(&mater, spt::matd<dim>{ first, second });
        #endif
    }

    for (std::size_t i = 0; i < init_poses.size(); ++i)
        automata.spawn_grain(&grains[i], automata.offset(init_poses[i]), cgr::nbh::nbhood_kind::euclid);

    
    //#define SHOWPIC
    progress_bar bar("Crystallization", automata.num_cells(), 70);
    while (!automata.stop_condition()) {
        //for (std::size_t i = 0; i < 5; ++i) {
        while (!automata.stop_condition()) {
            automata.iterate();
            bar.set_count(automata.num_crysted_cells());
        }
        
        #ifdef SHOWPIC
        std::ofstream ofile("automata-image-data.txt");
        ofile << "size " << size << std::endl;
        write_image_pixels(ofile, automata, false);
        ofile.close();
        std::system("python ./visualize.py");
        #endif // SHOWPIC
    }
    
    automata.smooth(4);
    #ifdef SHOWPIC
    std::ofstream ofile("automata-image-data.txt");
    ofile << "size " << size << std::endl;
    write_image_pixels(ofile, automata, false);
    ofile.close();
    std::system("python ./visualize.py");
    #endif // SHOWPIC

    #ifdef DIM3
    cgr::geo_from_automata simplegeo(&automata);
    simplegeo.make();

    std::cout << "planarity: " << simplegeo.planarity() << std::endl;
    std::cout << "before nonplanarity optimization:" << std::endl;
    std::cout << "  best nonplanarity/side_len: " << simplegeo.best_nonplanarity() / size << std::endl;
    std::cout << "  worst nonplanarity/side_len: " << simplegeo.worst_nonplanarity() / size << std::endl;
    std::cout << "  gmsh nonplanarity/side_len: " << simplegeo.gmsh_nonplanarity() / size << std::endl;
    simplegeo.optimize_nonplanarity();
    std::cout << "after nonplanarity optimization:" << std::endl;
    std::cout << "  gmsh nonplanarity/side_len: " << simplegeo.gmsh_nonplanarity() / size << std::endl;

    std::ofstream file("polycr.geo");
    simplegeo.write_geo(file);
    simplegeo.write_geo(std::cout);
    #endif // DIM3

    return 0;
}

int main() {
    seed = 2;
    for (std::size_t i = 0;; ++i) {
        std::cout << "iteration: " << i << std::endl;
        std::cout << "seed: " << seed << std::endl;
        try {
            inner_main();
        } catch (int _) {
            ++seed;
            continue;
        }
        break;
    }

    return 0;
}
