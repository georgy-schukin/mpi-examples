#include "neighbor.h"

namespace {
    std::vector<int> getDisplacements(bool has_prev, bool has_next) {
        if (has_prev && has_next) {
            return {-1, 0, 1};
        } else if (has_prev) {
            return {-1, 0};
        } else if (has_next) {
            return {0, 1};
        } else {
            return {0};
        }
    }
}

NeighborIndices::NeighborIndices() {
    generateIndices();
}

void NeighborIndices::generateIndices() {

}

const std::vector<Index3>& NeighborIndices::getNeighborIndices(int type) const {
    auto it = indices.find(type);
    if (it == indices.end()) {
        return indices.insert(IndexMap::value_type(type, buildIndices(type))).first->second;
    }
    return it->second;
}

std::vector<Index3> NeighborIndices::buildIndices(int type) {
    std::vector<Index3> neigh_indices;
    std::array<std::vector<int>, 3> displ;
    displ[0] = getDisplacements(type & PREV_X, type & NEXT_X);
    displ[1] = getDisplacements(type & PREV_Y, type & NEXT_Y);
    displ[2] = getDisplacements(type & PREV_Z, type & NEXT_Z);
    for (auto dx: displ[0])
    for (auto dy: displ[1])
    for (auto dz: displ[2]) {
        Index3 neigh_index {dx, dy, dz};
        if (dx == 0 && dy == 0 && dz == 0) {
            continue;
        }
        neigh_indices.push_back(neigh_index);
    }
    return neigh_indices;
}

int NeighborIndices::getType(bool px, bool nx, bool py, bool ny, bool pz, bool nz) {
    int type = 0;
    type |= (px ? PREV_X : 0);
    type |= (nx ? NEXT_X : 0);
    type |= (py ? PREV_Y : 0);
    type |= (ny ? NEXT_Y : 0);
    type |= (pz ? PREV_Z : 0);
    type |= (nz ? NEXT_Z : 0);
    return type;
}


