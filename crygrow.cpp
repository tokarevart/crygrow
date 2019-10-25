#include <iostream>
#include <fstream>
#include "simplest-automata.h"


constexpr std::size_t dim = 2;
using pos_t = spt::veci<dim>;
using automata_t = cgr::simplest_automata<dim>;
using cell_t = cgr::simplest_cell<dim>;
using crystallite_t = cgr::simplest_crystallite<dim>;
using material_t = cgr::simplest_material<dim>;
using nbhood_pos_t = cgr::nbhood_pos<dim>;

using pair_pos_cell = std::pair<std::vector<pos_t>, std::vector<cell_t>>;

std::vector<pos_t> make_poses_box(std::size_t size) {
    std::vector<pos_t> res;
    res.reserve(size * size);
    for (std::size_t i = 0; i < size; ++i)
        for (std::size_t j = 0; j < size; ++j)
            res.emplace_back(i, j);

    return res;
}


pair_pos_cell make_cells_box(std::size_t size, const cell_t& cell) {
    pair_pos_cell res;
    res.first = make_poses_box(size);
    res.second.assign(res.first.size(), cell);
    return res;
}


std::vector<pos_t> make_init_central_poses(std::int64_t ssize) {
    std::vector<pos_t> res;

    for (std::size_t i = 1; i <= 3; ++i)
        for (std::size_t j = 1; j <= 3; ++j)
            res.emplace_back(ssize * i / 4, ssize * j / 4);

    for (std::size_t i = 1; i <= 7; i += 2) {
        res.emplace_back((ssize * i) / 8, (ssize * 1) / 8);
        res.emplace_back(ssize * i / 8, ssize * 7 / 8);
    }
    res.emplace_back(ssize * 1 / 8, ssize * 3 / 8);
    res.emplace_back(ssize * 7 / 8, ssize * 3 / 8);
    res.emplace_back(ssize * 1 / 8, ssize * 5 / 8);
    res.emplace_back(ssize * 7 / 8, ssize * 5 / 8);

    return res;
}


int main() {
    std::size_t size = 400;
    std::int64_t ssize = size;
    std::size_t range = 7;
    cgr::nbhood_kind kind = cgr::nbhood_kind::euclid;
    automata_t automata({ 0, 0 }, { ssize, ssize }, range, kind);
    auto [default_poses, default_cells] = make_cells_box(size, cell_t());
    automata.set_cells(default_poses, default_cells);

    auto init_central_poses = make_init_central_poses(ssize);

    material_t mater;
    std::vector<crystallite_t> crysts(init_central_poses.size(), crystallite_t(&mater));

    std::vector<cell_t> init_central_cells;
    init_central_cells.reserve(crysts.size());
    for (auto& cryst : crysts)
        init_central_cells.emplace_back(1.0, &cryst);
    automata.set_cells(init_central_poses, init_central_cells);

    std::vector<nbhood_pos_t> init_nbhood_poses;
    init_nbhood_poses.reserve(init_central_poses.size());
    for (auto& pos : init_central_poses)
        init_nbhood_poses.emplace_back(cgr::make_nbhood_pos<2>(pos, kind, range));

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
            automata.iterate();

        for (std::size_t i = 0; i < automata.num_cells(); ++i) {
            auto pcell = automata.get_cell(i);
            if (pcell->crystallinity < 1.0 - automata_t::epsilon * (1.0 + pcell->crystallinity) ||
                pcell->crystallites.empty())
                continue;

            auto curpos = automata.pos(i);
            std::array<std::uint8_t, 3> color;
            if (pcell->crystallites.size() == 1) {
                color = { 0, 0, 0 };
            } else if (pcell->crystallites.size() == 2) {
                color = { 0, 0, 255 };
            } else {
                color = { 255, 0, 0 };
            }
            
            ofile 
                << curpos[0] << ' ' 
                << curpos[1] << ' '
                << static_cast<int>(color[0]) << ' ' 
                << static_cast<int>(color[1]) << ' ' 
                << static_cast<int>(color[2]) << std::endl;
        }
        ofile.close();
        std::system("python ./visualize.py");
    }

    return 0;
}
