#pragma once

#include "index.h"

#include <vector>
#include <map>

class NeighborIndices {
public:
    enum NeighborType {
        PREV_X = 0x1,
        NEXT_X = 0x2,
        PREV_Y = 0x4,
        NEXT_Y = 0x8,
        PREV_Z = 0x10,
        NEXT_Z = 0x20
    };

public:
    NeighborIndices();

    const std::vector<Index3>& getNeighborIndices(int type) const;

    static int getType(bool px, bool nx, bool py, bool ny, bool pz, bool nz);

private:
    void generateIndices();
    static std::vector<Index3> buildIndices(int type);

private:
    using IndexMap = std::map<int, std::vector<Index3>>;
    mutable IndexMap indices;
};
