#include <iostream>
#include <fstream>
#include "simplest-automata.h"


std::vector<cgr::simplest_cell<2>> set_cells_box(cgr::simplest_automata<2>& automata, std::size_t size) {
    std::vector<cgr::simplest_cell<2>> cells;
    automata.reserve(size * size);
    for (std::size_t i = 0; i < size; i++)
        for (std::size_t j = 0; j < size; j++) {
            cells.emplace_back();
            automata.set_cell({ i, j }, &cells.back());
        }
}


using pair_pos_cell = std::pair<std::vector<spt::veci<2>>, std::vector<cgr::simplest_cell<2>>>;

std::vector<spt::veci<2>> make_poses_box(std::size_t size) {
    std::vector<spt::veci<2>> res;
    res.reserve(size * size);
    for (std::size_t i = 0; i < size; i++)
        for (std::size_t j = 0; j < size; j++)
            res.emplace_back(i, j);

    return res;
}


pair_pos_cell make_cells_box(std::size_t size, const cgr::simplest_cell<2>& cell) {
    pair_pos_cell res;
    res.first = make_poses_box(size);
    res.second.assign(res.first.size(), cell);
    return res;
}


int main() {
    std::size_t size = 501;
    std::int64_t ssize2 = size / 2;
    cgr::simplest_automata<2> automata(10, cgr::nbhood_kind::euclid);
    auto [default_poses, default_cells] = make_cells_box(size, cgr::simplest_cell<2>());

    cgr::simplest_material<2> mater;
    cgr::simplest_crystallite<2> cryst(&mater);

    spt::veci<2> init_central_pos{ ssize2, ssize2 };
    cgr::simplest_cell<2> init_central_cell(1.0, &cryst);
    automata.set_cell(init_central_pos, &init_central_cell);

    auto init_neighbor_poses = cgr::make_nbhood_pos<2>(
        init_central_pos, cgr::nbhood_kind::euclid, 14);

    auto init_neighbor_cells = std::vector<cgr::simplest_cell<2>>(
        init_neighbor_poses.size(), init_central_cell);

    automata.set_cells(init_neighbor_poses, init_neighbor_cells);
    
    while (true) {
        std::ofstream ofile("automata-image-data.txt");
        ofile << "size " << size << std::endl;
        for (std::size_t i = 0; i < 40; i++)
            automata.iterate();

        for (std::size_t i = 0; i < automata.num_cells(); i++) {
            auto pcell = automata.get_cell(i);
            if (pcell->crystallinity == 0.0)
                continue;

            auto pos = automata.get_pos(i);
            auto brightness = static_cast<std::int64_t>(std::abs(1.0 - pcell->crystallinity) * 255.5);
            ofile << pos[0] << ' ' << pos[1] << ' '
                << brightness << ' ' << brightness << ' ' << brightness << std::endl;
        }
        ofile.close();
        std::system("python ./visualize.py");
    }

    return 0;
}
