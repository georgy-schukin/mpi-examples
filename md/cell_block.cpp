#include "cell_block.h"
#include "md.h"
#include "neighbor.h"

#include <map>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <sstream>

namespace {
    template <size_t Dim, bool AddShadow>
    Cell& getCell(CellBlock &block, const CellBlockParams &params, int fixed, int i1, int i2);

    template <>
    Cell& getCell<0, true>(CellBlock &block, const CellBlockParams &params, int fixed, int i1, int i2) {
        return block.cellAt(fixed + params.shadow_start[0], i1 + params.shadow_start[1], i2 + params.shadow_start[2]);
    }

    template <>
    Cell& getCell<1, true>(CellBlock &block, const CellBlockParams &params, int fixed, int i1, int i2) {
        return block.cellAt(i1 + params.shadow_start[0], fixed + params.shadow_start[1], i2 + params.shadow_start[2]);
    }

    template <>
    Cell& getCell<2, true>(CellBlock &block, const CellBlockParams &params, int fixed, int i1, int i2) {
        return block.cellAt(i1 + params.shadow_start[0], i2 + params.shadow_start[1], fixed + params.shadow_start[2]);
    }

    template <>
    Cell& getCell<0, false>(CellBlock &block, const CellBlockParams&, int fixed, int i1, int i2) {
        return block.cellAt(fixed, i1, i2);
    }

    template <>
    Cell& getCell<1, false>(CellBlock &block, const CellBlockParams&, int fixed, int i1, int i2) {
        return block.cellAt(i1, fixed, i2);
    }

    template <>
    Cell& getCell<2, false>(CellBlock &block, const CellBlockParams&, int fixed, int i1, int i2) {
        return block.cellAt(i1, i2, fixed);
    }

    template <size_t Dim, bool AddShadow>
    CellSlice2D getSlice(CellBlock &block, int fixed, int s1, int s2, int n1, int n2) {
        CellSlice2D slice(n1, n2);
        for (int i1 = 0; i1 < n1; i1++)
        for (int i2 = 0; i2 < n2; i2++) {
            slice.cellAt(i1, i2) = getCell<Dim, AddShadow>(block, block.getParams(), fixed, s1 + i1, s2 + i2);
        }
        return slice;
    }

    template <size_t Dim, bool AddShadow>
    void setSlice(CellBlock &block, int fixed, int s1, int s2, const CellSlice2D &slice) {
        for (int i1 = 0; i1 < slice.sizeX(); i1++)
        for (int i2 = 0; i2 < slice.sizeY(); i2++) {
             getCell<Dim, AddShadow>(block, block.getParams(), fixed, s1 + i1, s2 + i2) = slice.cellAt(i1, i2);
        }
    }
}

CellBlock::CellBlock(const CellBlockParams &params) :
    params(params) {
    initCells();
}

void CellBlock::setParams(const CellBlockParams &params) {
    this->params = params;
    initCells();
}

void CellBlock::initCells() {
    cells.resize(static_cast<size_t>(params.full_size[0] * params.full_size[1] * params.full_size[2]));
    forEachCellInd([this](Cell &cell, const Index3 &index) {
        cell.setIndex(Index3 {params.shift[0] + index[0] - params.shadow_start[0],
                              params.shift[1] + index[1] - params.shadow_start[1],
                              params.shift[2] + index[2] - params.shadow_start[2]});
        const auto neigh_type = NeighborIndices::getType(index[0] > 0, index[0] < params.full_size[0] - 1,
                                                         index[1] > 0, index[1] < params.full_size[1] - 1,
                                                         index[2] > 0, index[2] < params.full_size[2] - 1);
        cell.setNeighborIndicesType(neigh_type);
    });
}

Extent3 CellBlock::getExtent(const MeshParams &m_params) const {
    Extent3 extent;
    for (int i = 0; i < 3; i++) {
        extent.start(i) = m_params.origin[i] + params.shift[i] * m_params.step[i];
        extent.end(i) = extent.start(i) + params.size[i] * m_params.step[i];
    }
    return extent;
}

Cell& CellBlock::cellAt(int x, int y, int z) {
    return cells[static_cast<size_t>(cellIndex(x, y, z))];
}

const Cell& CellBlock::cellAt(int x, int y, int z) const {
    return cells[static_cast<size_t>(cellIndex(x, y, z))];
}

Cell& CellBlock::cellAt(const Index3 &index) {
    return cellAt(index[0], index[1], index[2]);
}

const Cell& CellBlock::cellAt(const Index3 &index) const {
    return cellAt(index[0], index[1], index[2]);
}

int CellBlock::cellIndex(int x, int y, int z) const {
    return x * params.full_size[1] * params.full_size[2] + y * params.full_size[2] + z;
}

CellSlice2D CellBlock::getSlice(int dimension, int x, int y, int z, int n1, int n2) {
    x = toIndex(x, 0);
    y = toIndex(y, 1);
    z = toIndex(z, 2);
    switch (dimension) {
        case 0: return ::getSlice<0, true>(*this, x, y, z, n1, n2);
        case 1: return ::getSlice<1, true>(*this, y, x, z, n1, n2);
        case 2: return ::getSlice<2, true>(*this, z, x, y, n1, n2);
    }
    return CellSlice2D();
}

void CellBlock::setSlice(int dimension, int x, int y, int z, const CellSlice2D &slice) {
    x = toIndex(x, 0);
    y = toIndex(y, 1);
    z = toIndex(z, 2);
    switch (dimension) {
        case 0: return ::setSlice<0, true>(*this, x, y, z, slice);
        case 1: return ::setSlice<1, true>(*this, y, x, z, slice);
        case 2: return ::setSlice<2, true>(*this, z, x, y, slice);
    }
}

void CellBlock::setSliceShadow(int dimension, int x, int y, int z, const CellSlice2D &slice) {
    x = toIndexFull(x, 0);
    y = toIndexFull(y, 1);
    z = toIndexFull(z, 2);
    switch (dimension) {
        case 0: return ::setSlice<0, false>(*this, x, y, z, slice);
        case 1: return ::setSlice<1, false>(*this, y, x, z, slice);
        case 2: return ::setSlice<2, false>(*this, z, x, y, slice);
    }
}

std::vector<Particle> CellBlock::extractAndMoveOutParticles(const MeshParams &mesh_params) {
    std::vector<Particle> out_particles;
    std::vector<Particle> in_particles;
    const auto extent = getExtent(mesh_params);
    forEachCell([mesh_params, extent, &in_particles, &out_particles](Cell &cell) {
        auto out_ps = cell.extractOutParticles(mesh_params);
        for (const auto &particle: out_ps) {
            if (particle.isIn(extent)) {
                in_particles.push_back(particle);
            } else {
                out_particles.push_back(particle);
            }
        }
    });
    addParticles(in_particles, mesh_params);
    return out_particles;
}

void CellBlock::addParticles(const std::vector<Particle> &particles, const MeshParams &mesh_params) {
    std::map<Index3, std::vector<Particle>> in_particles;
    const auto extent = getExtent(mesh_params);
    for (const auto &particle: particles) {
        const auto ix = int((particle.pos.x - extent.start(0)) / mesh_params.step[0]);
        const auto iy = int((particle.pos.y - extent.start(1)) / mesh_params.step[1]);
        const auto iz = int((particle.pos.z - extent.start(2)) / mesh_params.step[2]);
        if (ix < 0 || ix >= params.size[0] ||
            iy < 0 || iy >= params.size[1] ||
            iz < 0 || iz >= params.size[2]) {
            continue;
        }
        in_particles[Index3 {
                ix + params.shadow_start[0],
                iy + params.shadow_start[1],
                iz + params.shadow_start[2]}].push_back(particle);
    }
    for (const auto &p: in_particles) {
        cellAt(p.first).addParticles(p.second);
    }
}

int CellBlock::getNumOfParticles() const {
    int num_of_particles = 0;
    forEachCell([&num_of_particles](const Cell &cell) {
        num_of_particles += cell.getParticles().size();
    });
    return num_of_particles;
}

int CellBlock::toIndex(int index, int dim) const {
    return index >= 0 ? index : params.size[dim] + index; // negative indices count from end backwards
}

int CellBlock::toIndexFull(int index, int dim) const {
    return index >= 0 ? index : params.full_size[dim] + index; // negative indices count from end backwards
}

