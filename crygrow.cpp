#include <iostream>
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
    cgr::simplest_automata<2> automata;
    set_cells_box(automata, 301);
    std::vector<spt::vec<2>> grow_dirs;
    grow_dirs.emplace_back(-1.0, 1.0);
    grow_dirs.emplace_back(1.0, 1.0);
    grow_dirs.emplace_back(-1.0, -1.0);
    grow_dirs.emplace_back(1.0, -1.0);
    cgr::simplest_material<2> mater(grow_dirs);
    cgr::simplest_crystallite<2> cryst(&mater, spt::mat<2>{ 
        spt::vec<2>{ 1.0, 0.0 }, spt::vec<2>{ 0.0, 1.0 }});
    auto* init_cell = new cgr::simplest_cell<2>;
    init_cell->crystallinity = 1;
    init_cell->crystallites.assign(1, &cryst);
    automata.cell(spt::veci<2>{ 150, 150 }, init_cell);
    
    for (std::size_t i = 0; i < 10; i++)
        automata.iterate();

    return 0;
}
