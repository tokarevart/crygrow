#include <iostream>
#include <fstream>
#include "simplest-automata.h"


std::vector<cgr::simplest_cell<2>> set_cells_box(cgr::simplest_automata<2>& automata, std::size_t size) {
    std::vector<cgr::simplest_cell<2>> cells;
    automata.reserve(size * size);
    for (std::size_t i = 0; i < size; i++)
        for (std::size_t j = 0; j < size; j++) {
            cells.emplace_back();
            automata.set_cell({ 
                static_cast<std::int64_t>(i), 
                static_cast<std::int64_t>(j) }, &cells.back());
        }
}


int main() {
    std::size_t size = 501;
    std::int64_t ssize2 = size / 2;
    cgr::simplest_automata<2> automata(10, cgr::nbhood_kind::euclid);
    // TODO: simplify
    auto cells = set_cells_box(automata, size); 

    cgr::simplest_material<2> mater;
    cgr::simplest_crystallite<2> cryst(&mater);

    spt::veci<2> init_pos{ ssize2, ssize2 };

    cells.emplace_back(1.0, &cryst);
    automata.set_cell(init_pos, &cells.back());
    for (auto& pos : cgr::make_nbhood_pos<2>(init_pos, cgr::nbhood_kind::euclid, 14)) {
        cells.emplace_back(1.0, &cryst);
        automata.set_cell(pos, &cells.back());
    }
    
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
