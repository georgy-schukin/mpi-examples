#pragma once

#include "cell.h"
#include "cell_block.h"

#include "didal.h"

#include <string>
#include <functional>

struct MDParams {
    double cutoff_radius;
    double sigma;
    double epsilon;
    double delta_t;

    std::string output;
    int output_step = 1;

    bool periodic[3] = {false, false, false};

    enum ParticlesInitType {
        INIT_RANDOM = 0,
        INIT_RANDOM_ELLIPSE,
        INIT_COLLISION,
        INIT_EXPLOSION
    };

    ParticlesInitType particles_init = INIT_RANDOM;

    double init_speed = 1.0;
    double min_mass = 0.1;
    double max_mass = 1.0;

    bool load(const std::string &filename);
};

struct MeshParams {
    int mesh_size[3]; // num of cells in each dimension
    int num_blocks[3]; // num of cell blocks in each dimension
    double area_size[3]; // size of area
    double origin[3]; // origin - coordinate of (0, 0, 0) local point
    double step[3]; // cell size in each dimension
    bool periodic[3] = {false, false, false};

    Extent3 getExtent() const;
};

class NeighborIndices;

Vector3 calcForce(const Particle &p1, const Particle &p2, const MDParams &md_params);

void initForcesCell(Cell &cell);
void addForcesCell(Cell &cell, const MDParams &md_params);
void addForcesCell(Cell &dst_cell, const Cell &src_cell, const MDParams &md_params);
void updatePositionsCell(Cell &cell, double delta_t);
void updateVelocitiesCell(Cell &cell, double delta_t);

void calcForcesBlock(CellBlock &block, const MDParams &md_params, const NeighborIndices &neigh_indices);
void updatePositionsBlock(CellBlock &block, double delta_t);
void updateVelocitiesBlock(CellBlock &block, double delta_t);
