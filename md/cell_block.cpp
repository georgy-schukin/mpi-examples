#include "cell_block.h"
#include "md.h"
#include "neighbor.h"

#include <map>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <sstream>

using Range2D = ddl::IndexRangeND<2>;
using Range3D = ddl::IndexRangeND<3>;

namespace {
    Cell& cellAt(CellBlock::CellArray3 &cells, int x, int y, int z) {
        if (x < 0 || x >= static_cast<int>(cells.size(0)) ||
            y < 0 || y >= static_cast<int>(cells.size(1)) ||
            z < 0 || z >= static_cast<int>(cells.size(2))) {
            throw std::runtime_error("Access to non-existing cell (" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")");
        }
        return cells(x, y, z);
    }

    template <size_t Dim, bool AddShadow>
    Cell& getCell(CellBlock::CellArray3 &cells, const CellBlockParams &params, int fixed, int i1, int i2);

    template <>
    Cell& getCell<0, true>(CellBlock::CellArray3 &cells, const CellBlockParams &params, int fixed, int i1, int i2) {
        return cellAt(cells, fixed + params.shadow_start[0], i1 + params.shadow_start[1], i2 + params.shadow_start[2]);
    }

    template <>
    Cell& getCell<1, true>(CellBlock::CellArray3 &cells, const CellBlockParams &params, int fixed, int i1, int i2) {
        return cellAt(cells, i1 + params.shadow_start[0], fixed + params.shadow_start[1], i2 + params.shadow_start[2]);
    }

    template <>
    Cell& getCell<2, true>(CellBlock::CellArray3 &cells, const CellBlockParams &params, int fixed, int i1, int i2) {
        return cellAt(cells, i1 + params.shadow_start[0], i2 + params.shadow_start[1], fixed + params.shadow_start[2]);
    }

    template <>
    Cell& getCell<0, false>(CellBlock::CellArray3 &cells, const CellBlockParams&, int fixed, int i1, int i2) {
        return cellAt(cells, fixed, i1, i2);
    }

    template <>
    Cell& getCell<1, false>(CellBlock::CellArray3 &cells, const CellBlockParams&, int fixed, int i1, int i2) {
        return cellAt(cells, i1, fixed, i2);
    }

    template <>
    Cell& getCell<2, false>(CellBlock::CellArray3 &cells, const CellBlockParams&, int fixed, int i1, int i2) {
        return cellAt(cells, i1, i2, fixed);
    }

    template <size_t Dim, bool AddShadow>
    CellBlock::CellArray2 getSlice(CellBlock &block, int fixed, int s1, int s2, int n1, int n2) {
        try {
            CellBlock::CellArray2 slice(Range2D {static_cast<size_t>(n1), static_cast<size_t>(n2)});
            for (int i1 = 0; i1 < n1; i1++)
            for (int i2 = 0; i2 < n2; i2++) {
                slice(i1, i2) = getCell<Dim, AddShadow>(block.getCells(), block.getParams(), fixed, s1 + i1, s2 + i2);
            }
            return slice;
        }
        catch (std::bad_alloc&) {
            std::ostringstream out;
            const auto &cells = block.getCells();
            out << "Failed to alloc slice " << n1 << "x" << n2 <<
                   " from block " << cells.size(0) << "x" << cells.size(1) << "x" << cells.size(2) <<
                   " with " << block.getNumOfParticles() << " particles";
            throw std::runtime_error(out.str());
        }
    }

    template <size_t Dim, bool AddShadow>
    void setSlice(CellBlock &block, int fixed, int s1, int s2, const CellBlock::CellArray2 &slice) {
        for (int i1 = 0; i1 < static_cast<int>(slice.size(0)); i1++)
        for (int i2 = 0; i2 < static_cast<int>(slice.size(1)); i2++) {
             getCell<Dim, AddShadow>(block.getCells(), block.getParams(), fixed, s1 + i1, s2 + i2) = slice(i1, i2);
        }
    }

    template <size_t Dim, bool AddShadow>
    void setSlice(CellBlock &block, int fixed, int s1, int s2, CellBlock::CellArray2 &&slice) {
        for (int i1 = 0; i1 < static_cast<int>(slice.size(0)); i1++)
        for (int i2 = 0; i2 < static_cast<int>(slice.size(1)); i2++) {
             getCell<Dim, AddShadow>(block.getCells(), block.getParams(), fixed, s1 + i1, s2 + i2) = std::move(slice(i1, i2));
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
    cells = CellArray3(Range3D {static_cast<size_t>(params.full_size[0]),
                                static_cast<size_t>(params.full_size[1]),
                                static_cast<size_t>(params.full_size[2])});
    forEachCellInd([this](Cell &cell, const ddl::ArrayIndex<3> &index) {
        Index3 ind {static_cast<int>(index[0]),
                    static_cast<int>(index[1]),
                    static_cast<int>(index[2])};
        cell.setIndex(Index3 {params.shift[0] + ind[0] - params.shadow_start[0],
                              params.shift[1] + ind[1] - params.shadow_start[1],
                              params.shift[2] + ind[2] - params.shadow_start[2]});
        const auto neigh_type = NeighborIndices::getType(ind[0] > 0, ind[0] < params.full_size[0] - 1,
                                                         ind[1] > 0, ind[1] < params.full_size[1] - 1,
                                                         ind[2] > 0, ind[2] < params.full_size[2] - 1);
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

CellBlock::CellArray2 CellBlock::getSlice(int dimension, int x, int y, int z, int n1, int n2) {
    x = toIndex(x, 0);
    y = toIndex(y, 1);
    z = toIndex(z, 2);
    switch (dimension) {
        case 0: return ::getSlice<0, true>(*this, x, y, z, n1, n2);
        case 1: return ::getSlice<1, true>(*this, y, x, z, n1, n2);
        case 2: return ::getSlice<2, true>(*this, z, x, y, n1, n2);
    }
    return CellArray2();
}

void CellBlock::setSlice(int dimension, int x, int y, int z, const CellArray2 &slice) {
    x = toIndex(x, 0);
    y = toIndex(y, 1);
    z = toIndex(z, 2);
    switch (dimension) {
        case 0: return ::setSlice<0, true>(*this, x, y, z, slice);
        case 1: return ::setSlice<1, true>(*this, y, x, z, slice);
        case 2: return ::setSlice<2, true>(*this, z, x, y, slice);
    }
}

void CellBlock::setSlice(int dimension, int x, int y, int z, CellArray2 &&slice) {
    x = toIndex(x, 0);
    y = toIndex(y, 1);
    z = toIndex(z, 2);
    switch (dimension) {
        case 0: return ::setSlice<0, true>(*this, x, y, z, std::move(slice));
        case 1: return ::setSlice<1, true>(*this, y, x, z, std::move(slice));
        case 2: return ::setSlice<2, true>(*this, z, x, y, std::move(slice));
    }
}

void CellBlock::setSliceShadow(int dimension, int x, int y, int z, const CellArray2 &slice) {
    x = toIndexFull(x, 0);
    y = toIndexFull(y, 1);
    z = toIndexFull(z, 2);
    switch (dimension) {
        case 0: return ::setSlice<0, false>(*this, x, y, z, slice);
        case 1: return ::setSlice<1, false>(*this, y, x, z, slice);
        case 2: return ::setSlice<2, false>(*this, z, x, y, slice);
    }
}

void CellBlock::setSliceShadow(int dimension, int x, int y, int z, CellArray2 &&slice) {
    x = toIndexFull(x, 0);
    y = toIndexFull(y, 1);
    z = toIndexFull(z, 2);
    switch (dimension) {
        case 0: return ::setSlice<0, false>(*this, x, y, z, std::move(slice));
        case 1: return ::setSlice<1, false>(*this, y, x, z, std::move(slice));
        case 2: return ::setSlice<2, false>(*this, z, x, y, std::move(slice));
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
    std::map<ddl::ArrayIndex<3>, std::vector<Particle>> in_particles;
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
        in_particles[ddl::ArrayIndex<3> {
                static_cast<size_t>(ix + params.shadow_start[0]),
                static_cast<size_t>(iy + params.shadow_start[1]),
                static_cast<size_t>(iz + params.shadow_start[2])}].push_back(particle);
    }
    for (const auto &p: in_particles) {
        cells[p.first].addParticles(p.second);
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

