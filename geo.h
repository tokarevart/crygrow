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
    T* entity;
    orientation orient;

    void reverse() {
        orient = orient == orientation::forward ?
            orientation::reverse : orientation::forward;
    }

    auto begin() {
        return orient == orientation::forward ?
            entity->begin() : entity->rbegin();
    }
    auto end() {
        return orient == orientation::forward ?
            entity->end() : entity->rend();
    }
    auto rbegin() {
        return orient == orientation::forward ?
            entity->rbegin() : entity->begin();
    }
    auto rend() {
        return orient == orientation::forward ?
            entity->rend() : entity->end();
    }

    auto front_ptr() {
        return entity->front_ptr();
    }
    auto back_ptr() {
        return entity->back_ptr();
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
    std::array<point*, 2> points;

    auto begin() { 
        return points.begin(); 
    }
    auto end() { 
        return points.end(); 
    }
    auto rbegin() { 
        return points.rbegin(); 
    }
    auto rend() { 
        return points.rend(); 
    }

    auto front_ptr() {
        return points.front();
    }
    auto back_ptr() {
        return points.back();
    }
};

struct plane_surface {
    std::size_t tag;
    std::vector<std::unique_ptr<oriented<line>>> lines;

    auto begin() { 
        return lines.begin(); 
    }
    auto end() { 
        return lines.end(); 
    }
    auto rbegin() { 
        return lines.rbegin(); 
    }
    auto rend() { 
        return lines.rend(); 
    }

    auto front_ptr() { 
        return lines.front().get(); 
    }
    auto back_ptr() {
        return lines.back().get();
    }
};

struct surface {
    std::size_t tag;
    std::vector<std::unique_ptr<oriented<line>>> lines;

    auto begin() { 
        return lines.begin(); 
    }
    auto end() { 
        return lines.end(); 
    }
    auto rbegin() { 
        return lines.rbegin(); 
    }
    auto rend() { 
        return lines.rend(); 
    }

    auto front_ptr() {
        return lines.front().get();
    }
    auto back_ptr() {
        return lines.back().get();
    }
};

struct volume {
    std::size_t tag;
    std::vector<std::unique_ptr<oriented<surface>>> surfaces;

    auto begin() { 
        return surfaces.begin(); 
    }
    auto end() { 
        return surfaces.end(); 
    }
    auto rbegin() { 
        return surfaces.rbegin(); 
    }
    auto rend() { 
        return surfaces.rend(); 
    }

    auto front_ptr() {
        return surfaces.front().get();
    }
    auto back_ptr() {
        return surfaces.back().get();
    }
};

}
