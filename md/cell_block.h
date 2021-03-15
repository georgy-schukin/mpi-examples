#pragma once

#include "cell.h"
#include "cell_slice_2d.h"
#include "index.h"

#include <vector>

struct CellBlockParams {
    int size[3]; // num of cells in each dimension
    int shift[3]; // shift from start in cells in each dimension
    int shadow_start[3];
    int shadow_end[3];
    int full_size[3]; // total num of cells in each dimenison
    int neighbor_indices_type = 0;
};

struct MeshParams;

class CellBlock {
public:
    CellBlock() {}
    CellBlock(const CellBlockParams &params);

    std::vector<Cell>& getCells() {
        return cells;
    }

    const std::vector<Cell>& getCells() const {
        return cells;
    }

    void setParams(const CellBlockParams &params);

    const CellBlockParams& getParams() const {
        return params;
    }

    Extent3 getExtent(const MeshParams &m_params) const;

    Cell& cellAt(int x, int y, int z);
    const Cell& cellAt(int x, int y, int z) const;

    Cell& cellAt(const Index3 &index);
    const Cell& cellAt(const Index3 &index) const;

    int cellIndex(int x, int y, int z) const;

    CellSlice2D getSlice(int dimension, int x, int y, int z, int n1, int n2);

    void setSlice(int dimension, int x, int y, int z, const CellSlice2D &slice);
    void setSliceShadow(int dimension, int x, int y, int z, const CellSlice2D &slice);

    std::vector<Particle> extractAndMoveOutParticles(const MeshParams &mesh_params);
    void addParticles(const std::vector<Particle> &particles, const MeshParams &mesh_params);

    template <typename Proc, typename... Args>
    void forEachCell(Proc proc, Args... args) {
        const auto p = getParams();
        for (int i = p.shadow_start[0]; i < p.full_size[0] - p.shadow_end[0]; i++)
        for (int j = p.shadow_start[1]; j < p.full_size[1] - p.shadow_end[1]; j++)
        for (int k = p.shadow_start[2]; k < p.full_size[2] - p.shadow_end[2]; k++) {
            proc(cellAt(i, j, k), args...);
        }
    }

    template <typename Proc, typename... Args>
    void forEachCell(Proc proc, Args... args) const {
        const auto p = getParams();
        for (int i = p.shadow_start[0]; i < p.full_size[0] - p.shadow_end[0]; i++)
        for (int j = p.shadow_start[1]; j < p.full_size[1] - p.shadow_end[1]; j++)
        for (int k = p.shadow_start[2]; k < p.full_size[2] - p.shadow_end[2]; k++) {
            proc(cellAt(i, j, k), args...);
        }
    }

    template <typename Proc, typename... Args>
    void forEachCellInd(Proc proc, Args... args) {
        const auto p = getParams();
        for (int i = p.shadow_start[0]; i < p.full_size[0] - p.shadow_end[0]; i++)
        for (int j = p.shadow_start[1]; j < p.full_size[1] - p.shadow_end[1]; j++)
        for (int k = p.shadow_start[2]; k < p.full_size[2] - p.shadow_end[2]; k++) {
            const auto index = Index3 {i, j, k};
            proc(cellAt(i, j, k), index, args...);
        }
    }

    template <typename Proc, typename... Args>
    void forEachCellInd(Proc proc, Args... args) const {
        const auto p = getParams();
        for (int i = p.shadow_start[0]; i < p.full_size[0] - p.shadow_end[0]; i++)
        for (int j = p.shadow_start[1]; j < p.full_size[1] - p.shadow_end[1]; j++)
        for (int k = p.shadow_start[2]; k < p.full_size[2] - p.shadow_end[2]; k++) {
            const auto index = Index3 {i, j, k};
            proc(cellAt(i, j, k), index, args...);
        }
    }

    template <typename Proc, typename... Args>
    void forEachParticle(Proc proc, Args... args) {
        forEachCell([proc, args...](Cell &cell) {
            cell.forEachParticle(proc, args...);
        });
    }

    template <typename Proc, typename... Args>
    void forEachParticle(Proc proc, Args... args) const {
        forEachCell([proc, args...](const Cell &cell) {
            cell.forEachParticle(proc, args...);
        });
    }

    int getNumOfParticles() const;

private:
    int toIndex(int index, int dim) const;
    int toIndexFull(int index, int dim) const;
    void initCells();

private:
    CellBlockParams params;
    std::vector<Cell> cells;
};
