#include <iostream>
#include <fstream>
#include "simplest-automata.h"


using pair_pos_cell = std::pair<std::vector<spt::veci<2>>, std::vector<cgr::simplest_cell<2>>>;

std::vector<spt::veci<2>> make_poses_box(std::size_t size) {
    std::vector<spt::veci<2>> res;
    res.reserve(size * size);
    for (std::size_t i = 0; i < size; ++i)
        for (std::size_t j = 0; j < size; ++j)
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
    std::size_t size = 401;
    std::int64_t ssize = size;
    std::size_t range = 7;
    cgr::nbhood_kind kind = cgr::nbhood_kind::euclid;
    cgr::simplest_automata<2> automata({ 0, 0 }, { ssize, ssize }, range, kind);
    auto [default_poses, default_cells] = make_cells_box(size, cgr::simplest_cell<2>());
    automata.set_cells(default_poses, default_cells);

    std::size_t num_inits = 21;
    cgr::simplest_material<2> mater;
    std::vector<cgr::simplest_crystallite<2>> crysts;
    crysts.reserve(num_inits);
    for (std::size_t i = 0; i < num_inits; ++i)
        crysts.emplace_back(&mater);

    std::vector<spt::veci<2>> init_central_poses;
    init_central_poses.reserve(num_inits);

    for (std::size_t i = 1; i <= 3; ++i)
        for (std::size_t j = 1; j <= 3; ++j)
            init_central_poses.emplace_back(ssize * i / 4, ssize * j / 4);

    for (std::size_t i = 1; i <= 7; i += 2) {
        init_central_poses.emplace_back((ssize * i) / 8, (ssize * 1) / 8);
        init_central_poses.emplace_back(ssize * i / 8, ssize * 7 / 8);
    }
    init_central_poses.emplace_back(ssize * 1 / 8, ssize * 3 / 8);
    init_central_poses.emplace_back(ssize * 7 / 8, ssize * 3 / 8);
    init_central_poses.emplace_back(ssize * 1 / 8, ssize * 5 / 8);
    init_central_poses.emplace_back(ssize * 7 / 8, ssize * 5 / 8);

    std::vector<cgr::simplest_cell<2>> init_central_cells;
    init_central_cells.reserve(num_inits);
    for (std::size_t i = 0; i < num_inits; ++i)
        init_central_cells.emplace_back(1.0, &crysts[i]);
    automata.set_cells(init_central_poses, init_central_cells);

    std::vector<cgr::nbhood_pos<2>> init_nbhood_poses;
    init_nbhood_poses.reserve(num_inits);
    for (std::size_t i = 0; i < num_inits; ++i)
        init_nbhood_poses.emplace_back(std::move(cgr::make_nbhood_pos<2>(
            init_central_poses[i], kind, range)));

    std::vector<std::vector<cgr::simplest_cell<2>>> init_nbhood_cells;
    init_nbhood_cells.reserve(num_inits);
    for (std::size_t i = 0; i < num_inits; ++i)
        init_nbhood_cells.emplace_back(init_nbhood_poses[i].size(), init_central_cells[i]);

    for (std::size_t i = 0; i < num_inits; ++i)
        automata.set_cells(init_nbhood_poses[i], init_nbhood_cells[i]);
    
    while (!automata.stop_condition()) {
        std::ofstream ofile("automata-image-data.txt");
        ofile << "size " << size << std::endl;
        for (std::size_t i = 0; i < 130; ++i)
            automata.iterate();

        for (std::size_t i = 0; i < automata.num_cells(); ++i) {
            auto pcell = automata.get_cell(i);
            if (pcell->crystallinity < 1.0 - cgr::simplest_automata<2>::epsilon ||
                pcell->crystallites.size() > 1)
                continue;

            auto curpos = automata.pos(i);
            auto brightness = static_cast<std::int64_t>(std::abs(1.0 - pcell->crystallinity) * 255.5);
            ofile << curpos[0] << ' ' << curpos[1] << ' '
                << brightness << ' ' << brightness << ' ' << brightness << std::endl;
        }
        ofile.close();
        std::system("python ./visualize.py");
    }

    return 0;
}
