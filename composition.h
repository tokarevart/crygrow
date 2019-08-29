#pragma once
#include "automata-base.h"

namespace cgr {

template <typename Automata>
class compose : automata_base<Automata::dim, Automata::cell_type> {
protected:
    Automata lower_level;
};

} // namespace cgr
