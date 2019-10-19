#include <iostream>
#include <fstream>
#include "simplest-automata.h"


void set_cells_box(cgr::simplest_automata<2>& automata, std::size_t size) {
    automata.reserve(size * size);
    for (std::size_t i = 0; i < size; i++)
        for (std::size_t j = 0; j < size; j++)
            automata.cell({ static_cast<std::int64_t>(i), 
                          static_cast<std::int64_t>(j) }, 
                          new cgr::simplest_cell<2>);
}


int main() {
    std::size_t size = 501;
    std::int64_t ssize2 = size / 2;
    cgr::simplest_automata<2> automata(10, cgr::nbhood_kind::euclid);
    set_cells_box(automata, size);
    std::vector<spt::vec<2>> grow_dirs;
    grow_dirs.emplace_back(-1.0, 1.0);
    grow_dirs.emplace_back(1.0, 1.0);
    grow_dirs.emplace_back(-1.0, -1.0);
    grow_dirs.emplace_back(1.0, -1.0);
    cgr::simplest_material<2> mater(grow_dirs);
    cgr::simplest_crystallite<2> cryst(&mater, spt::mat<2>{ 
        spt::vec<2>{ 1.0, 0.0 }, spt::vec<2>{ 0.0, 1.0 }});
    spt::veci<2> init_pos{ ssize2, ssize2 };
    automata.cell(init_pos, new cgr::simplest_cell<2>(1.0, &cryst));
    for (auto& pos : cgr::neighbors_pos(init_pos, cgr::nbhood_kind::euclid, 14))
        automata.cell(pos, new cgr::simplest_cell<2>(1.0, &cryst));
    
    while (true) {
        std::ofstream ofile("automata-image-data.txt");
        ofile << "size " << size << std::endl;
        for (std::size_t i = 0; i < 40; i++)
            automata.iterate();

        for (auto pcell : automata) {
            if (pcell->crystallinity == 0.0)
                continue;

            auto pos = automata.pos(pcell);
            auto brightness = static_cast<std::int64_t>(std::abs(1.0 - pcell->crystallinity) * 255.5);
            ofile << pos[0] << ' ' << pos[1] << ' '
                << brightness << ' ' << brightness << ' ' << brightness << std::endl;
        }
        ofile.close();
        std::system("python ./vis.py");
    }

    return 0;
}
