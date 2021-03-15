#pragma once

#include "cell.h"

#include <vector>

class CellSlice2D {
public:
    CellSlice2D() {}
    CellSlice2D(int sz_x, int sz_y) :
        size_x(sz_x), size_y(sz_y) {
        cells.resize(static_cast<size_t>(size_x * size_y));
    }

    int sizeX() const {
        return size_x;
    }

    int sizeY() const {
        return size_y;
    }

    Cell& cellAt(int x, int y) {
        return cells[x * size_y + y];
    }

    const Cell& cellAt(int x, int y) const {
        return cells[x * size_y + y];
    }

private:
    std::vector<Cell> cells;
    int size_x, size_y;
};
