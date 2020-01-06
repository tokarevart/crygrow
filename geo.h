#pragma once
#include <cstddef>
#include <array>
#include <vector>
#include <memory>
#include "vec.h"
#include "sptops.h"
#include "sptalgs.h"

namespace geo {

using real_type = double;

enum class orientation {
    forward,
    reverse
};

template <typename T>
struct oriented {
    const T* entity;
    orientation orient;

    void reverse() {
        orient = orient == orientation::forward ?
            orientation::reverse : orientation::forward;
    }

    oriented(const T* entity, orientation orient = orientation::forward)
        : entity(entity), orient(orient) {}
};

struct geometry {
    std::vector<std::unique_ptr<volume>>  volumes;
    std::vector<std::unique_ptr<surface>> surfaces;
    std::vector<std::unique_ptr<line>>    lines;
    std::vector<std::unique_ptr<point>>   points;
};

struct point {
    std::size_t tag;
    spt::vec3<real_type> x;
};

struct line {
    std::size_t tag;
    std::pair<point*, point*> points;
};

struct plane_surface {
    std::size_t tag;
    std::vector<oriented<line>> lines;
};

struct surface {
    std::size_t tag;
    std::vector<oriented<line>> lines;
};

struct volume {
    std::size_t tag;
    std::vector<oriented<surface>> surfaces;
};

}
