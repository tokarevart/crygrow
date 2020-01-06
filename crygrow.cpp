#include <iostream>
#include <fstream>
#include <random>
#include "sptalgs.h"
#include "simple-automata.h"
#include "simple-geometry.h"

#define DIM3

#ifdef DIM3
constexpr std::size_t dim = 3;
#else
constexpr std::size_t dim = 2;
#endif
std::size_t seed = 2;
constexpr auto kind = cgr::nbhood_kind::euclid;
using pos_t = spt::veci<dim>;
using automata_t = cgr::simple_automata<dim, kind>;
using cell_t = cgr::simple_cell<dim>;
using crystallite_t = cgr::simple_grain<dim>;
using material_t = cgr::simple_material<dim>;
using nbhood_pos_t = cgr::nbhood_pos<dim>;

using pair_pos_cell = std::pair<std::vector<pos_t>, std::vector<cell_t>>;


std::vector<pos_t> make_poses_box(std::size_t size) {
    std::vector<pos_t> res;
    std::size_t ressize = 1;
    for (std::size_t i = 0; i < dim; ++i)
        ressize *= size;
    res.reserve(ressize);
    for (std::size_t i = 0; i < ressize; ++i)
        res.emplace_back(cgr::upos(i, spt::vecu<dim>::filled_with(size)));

    //for (std::size_t i = 0; i < size; ++i)
    //    for (std::size_t j = 0; j < size; ++j)
    //        res.emplace_back(i, j);

    return res;
}

pair_pos_cell make_cells_box(std::size_t size, const cell_t& cell) {
    pair_pos_cell res;
    res.first = make_poses_box(size);
    res.second.assign(res.first.size(), cell);
    return res;
}

std::vector<pos_t> make_central_pos(std::size_t size) {
    std::int64_t ssize = size;
    std::vector<pos_t> res;
    res.emplace_back(pos_t::filled_with(ssize * 1 / 2));
    return res;
}

//std::vector<pos_t> make_central_poses(std::size_t size) {
//    std::int64_t ssize = size;
//
//    std::vector<pos_t> res;
//
//    for (std::size_t i = 1; i <= 3; ++i)
//        for (std::size_t j = 1; j <= 3; ++j)
//            res.emplace_back(ssize * i / 4, ssize * j / 4);
//
//    for (std::size_t i = 1; i <= 7; i += 2) {
//        res.emplace_back(ssize * i / 8, ssize * 1 / 8);
//        res.emplace_back(ssize * i / 8, ssize * 7 / 8);
//    }
//    res.emplace_back(ssize * 1 / 8, ssize * 3 / 8);
//    res.emplace_back(ssize * 7 / 8, ssize * 3 / 8);
//    res.emplace_back(ssize * 1 / 8, ssize * 5 / 8);
//    res.emplace_back(ssize * 7 / 8, ssize * 5 / 8);
//
//    return res;
//}

std::uint64_t min_distance2(pos_t pos, std::vector<pos_t> others) {
    std::uint64_t res = std::numeric_limits<std::uint64_t>::max();
    for (auto& other : others) {
        std::uint64_t dist = (other - pos).magnitude2();
        if (dist < res)
            res = dist;
    }
    return res;
}

std::vector<pos_t> make_random_central_poses(std::size_t size, std::size_t num, std::uint64_t min_dist2 = 0) {
    std::vector<pos_t> res;
    std::mt19937_64 gen(seed);
    std::uniform_int_distribution<std::size_t> dis(0, size - 1);
    for (std::size_t i = 0; i < num;) {
        pos_t curpos;
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

#include "iteration.h"
int main_test() {
    itr::iteration<int> it(0, 7);
    for (auto i = it.begin(); !it.has_ended(); i = it.next())
        std::cout << i << ' ';
    std::cout << std::endl;
    it.init(0, 7, itr::direction::reverse);
    for (auto i = it.begin(); !it.has_ended(); i = it.next())
        std::cout << i << ' ';

    return 0;
}

int main() {
    std::size_t size = 100;
    std::size_t range = 3;
    automata_t automata(size, range);
    auto [default_poses, default_cells] = make_cells_box(size, cell_t());
    automata.set_cells(default_poses, default_cells);

    //auto init_central_poses = make_central_pos(size);
    auto init_central_poses = make_random_central_poses(size, 5, (range * 4) * (range * 4));

    material_t mater(cgr::material_property::isotropic, { 
        #ifdef DIM3
        spt::vecd<dim>{ 1.0, 0.0, 0.0 }.normalize(),
        spt::vecd<dim>{ 0.0, 1.0, 0.0 }.normalize(),
        spt::vecd<dim>{ 0.0, 0.0, 1.0 }.normalize() });
        #else
        spt::vecd<dim>{ 4.0, 1.0 }.normalize(), 
        spt::vecd<dim>{ -1.0, 4.0 }.normalize() });
        #endif
    //std::vector<crystallite_t> crysts(init_central_poses.size(), crystallite_t(&mater));
    std::vector<crystallite_t> crysts;
    crysts.reserve(init_central_poses.size());
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> dis(-1.0, std::nextafter(1.0, 2.0));
    for (std::size_t i = 0; i < init_central_poses.size(); ++i) {
        #ifdef DIM3
        const double pi = 3.14159265359;
        auto rot = spt::rotation(spt::vecd<dim>{ dis(gen), dis(gen), dis(gen) }.normalize(), std::abs(dis(gen)) * pi);
        crysts.emplace_back(&mater, rot);
        #else
        auto first = spt::vecd<dim>{ dis(gen), dis(gen) }.normalize();
        spt::vecd<dim> second{ -first[1], first[0] };
        crysts.emplace_back(&mater, spt::matd<dim>{ first, second });
        #endif
    }

    std::vector<cell_t> init_central_cells;
    init_central_cells.reserve(crysts.size());
    for (auto& cryst : crysts)
        init_central_cells.emplace_back(1.0, &cryst);
    automata.set_cells(init_central_poses, init_central_cells);

    std::vector<nbhood_pos_t> init_nbhood_poses;
    init_nbhood_poses.reserve(init_central_poses.size());
    for (auto& pos : init_central_poses)
        init_nbhood_poses.emplace_back(cgr::make_nbhood_pos<kind, dim>(pos, range));

    std::vector<std::vector<cell_t>> init_nbhood_cells;
    init_nbhood_cells.reserve(init_nbhood_poses.size());
    for (std::size_t i = 0; i < init_nbhood_poses.size(); ++i)
        init_nbhood_cells.emplace_back(init_nbhood_poses[i].size(), init_central_cells[i]);

    for (std::size_t i = 0; i < init_nbhood_poses.size(); ++i)
        automata.set_cells(init_nbhood_poses[i], init_nbhood_cells[i]);
    
    while (!automata.stop_condition()) {
        std::ofstream ofile("automata-image-data.txt");
        ofile << "size " << size << std::endl;
        for (std::size_t i = 0; i < 100; ++i)
        while (!automata.stop_condition())
            automata.iterate();
        
        for (std::size_t i = 0; i < size * size; ++i) {
            auto curpos = automata.upos(i);
            #ifdef DIM3
            //curpos[2] = static_cast<std::int64_t>(size) / 2;
            //curpos[2] = size - 1;
            curpos[2] = 0;
            #endif
            auto pcell = automata.get_cell(curpos);

            std::array<std::uint8_t, 3> color;
            constexpr bool blackwhite = false;
            if constexpr (!blackwhite) {
                if (pcell->crystallinity < 1.0 - automata_t::epsilon * (1.0 + pcell->crystallinity) &&
                    !pcell->grains.empty()) {
                    color = { 0, 255, 0 };
                } else if (pcell->crystallinity < 1.0 - automata_t::epsilon * (1.0 + pcell->crystallinity) ||
                           pcell->grains.empty()) {
                    color = { 255, 255, 255 };
                } else if (pcell->grains.size() == 1) {
                    color = { 0, 0, 0 };
                } else if (pcell->grains.size() == 2) {
                    color = { 0, 0, 255 };
                } else {
                    color = { 255, 0, 0 };
                }
            } else {
                if (pcell->grains.size() == 1 &&
                    std::abs(pcell->crystallinity - 1.0) <= automata_t::epsilon * (1.0 + pcell->crystallinity)) {
                    color = { 0, 0, 0 };
                } else {
                    color = { 255, 255, 255 };
                }
            }
            
            
            ofile 
                << curpos[0] << ' ' 
                << curpos[1] << ' '
                << static_cast<int>(color[0]) << ' ' 
                << static_cast<int>(color[1]) << ' ' 
                << static_cast<int>(color[2]) << std::endl;
        }
        ofile.close();

        std::size_t num_cells_gtdimpl1 = 0;
        for (std::size_t i = 0; i < size * size; ++i) {
            auto pcell = automata.get_cell(i);
            if (pcell->grains.size() > dim + 1)
                ++num_cells_gtdimpl1;
        }
        std::cout << "number of cells with grains number >= dim + 1: " << num_cells_gtdimpl1 << std::endl;

        std::system("python ./visualize.py");
    }

    cgr::simple_geo_tool simplegeo(&automata);
    simplegeo.make_geometry();
    std::ofstream file("kek.geo");
    simplegeo.write_geo(file);
    simplegeo.write_geo(std::cout);

    return 0;
}
