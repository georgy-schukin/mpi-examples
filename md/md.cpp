#include "md.h"
#include "neighbor.h"

#include <random>
#include <fstream>
#include <cmath>
#include <limits>
#include <map>

bool MDParams::load(const std::string &filename) {
    static const std::map<std::string, ParticlesInitType> init_types = {
        {"random", INIT_RANDOM},
        {"random_ellipse", INIT_RANDOM_ELLIPSE},
        {"collision", INIT_COLLISION},
        {"explosion", INIT_EXPLOSION}
    };
    std::ifstream in(filename.c_str());
    if (!in.is_open()) {
        return false;
    }
    while (!in.eof()) {
        std::string tag;
        in >> tag;
        if (tag == "cutoff_r") {
            in >> cutoff_radius;
        } else if (tag == "sigma") {
            in >> sigma;
        } else if (tag == "epsilon") {
            in >> epsilon;
        } else if (tag == "delta_t") {
            in >> delta_t;
        } else if (tag == "output") {
            in >> output;
        } else if (tag == "output_step") {
            in >> output_step;
        } else if (tag == "periodic") {
            in >> periodic[0] >> periodic[1] >> periodic[2];
        } else if (tag == "init") {
            std::string init_type;
            in >> init_type;
            auto it = init_types.find(init_type);
            if (it != init_types.end()) {
                particles_init = it->second;
            }
        } else if (tag == "speed") {
            in >> init_speed;
        } else if (tag == "min_mass") {
            in >> min_mass;
        } else if (tag == "max_mass") {
            in >> max_mass;
        }
    }
    return true;
}

Extent3 MeshParams::getExtent() const {
    Extent3 extent;
    for (int i = 0; i < 3; i++) {
        extent.start(i) = origin[i];
        extent.end(i) = extent.start(i) + mesh_size[i] * step[i];
    }
    return extent;
}

Vector3 calcForce(const Particle &p1, const Particle &p2, const MDParams &md_params) {
    // Lennard-Jones force.
    Vector3 r(p2.pos.x - p1.pos.x,
              p2.pos.y - p1.pos.y,
              p2.pos.z - p1.pos.z);
    const double dist_squared = r.x * r.x + r.y * r.y + r.z * r.z + 1e-12;
    if (dist_squared > md_params.cutoff_radius * md_params.cutoff_radius) {
        return Vector3();
    }
    const double s = md_params.sigma * md_params.sigma / dist_squared;
    const double s3 = s * s * s; // (sigma / r) ^ 6
    const double f = 24.0 * md_params.epsilon * s3 * (1.0 - 2.0 * s3) / dist_squared;
    return Vector3(r.x * f, r.y * f, r.z * f);
}

void initForcesCell(Cell &cell) {
    for (auto &particle: cell.getParticles()) {
        particle.force.null();
    }
}

void addForcesCell(Cell &cell, const MDParams &md_params) {
    auto &particles = cell.getParticles();
    for (size_t i = 0; i < particles.size(); i++) {
        for (size_t j = i + 1; j < particles.size(); j++) {
            const auto force = calcForce(particles[i], particles[j], md_params);
            particles[i].force.add(force);
            particles[j].force.sub(force);
        }
    }
}

void addForcesCell(Cell &dst_cell, const Cell &src_cell, const MDParams &md_params) {
    for (auto &dst_particle: dst_cell.getParticles()) {
        for (const auto &src_particle: src_cell.getParticles()) {
            const auto force = calcForce(dst_particle, src_particle, md_params);
            dst_particle.force.add(force);
        }
    }
}

void updatePositionsCell(Cell &cell, double delta_t) {
    for (auto &particle: cell.getParticles()) {
        const auto dd = delta_t * delta_t / (2 * particle.mass);
        particle.pos.x += delta_t * particle.velocity.x + particle.force.x * dd;
        particle.pos.y += delta_t * particle.velocity.y + particle.force.y * dd;
        particle.pos.z += delta_t * particle.velocity.z + particle.force.z * dd;
        particle.force_prev = particle.force;
    }
}

void updateVelocitiesCell(Cell &cell, double delta_t) {
    for (auto &particle: cell.getParticles()) {
        const auto dd = delta_t / (2 * particle.mass);
        particle.velocity.x += (particle.force.x + particle.force_prev.x) * dd;
        particle.velocity.y += (particle.force.y + particle.force_prev.y) * dd;
        particle.velocity.z += (particle.force.z + particle.force_prev.z) * dd;
    }
}

void calcForcesBlock(CellBlock &block, const MDParams &md_params, const NeighborIndices &neigh_indices) {
    block.forEachCellInd([&block, &md_params, &neigh_indices](Cell &cell, const ddl::ArrayIndex<3> &index) {
        initForcesCell(cell);
        addForcesCell(cell, md_params);
        for (const auto &ni: neigh_indices.getNeighborIndices(cell.getNeighborIndicesType())) {
            ddl::ArrayIndex<3> neigh_index {index[0] + ni[0], index[1] + ni[1], index[2] + ni[2]};
            addForcesCell(cell, block.getCell(neigh_index), md_params);
        }
    });
}

void updatePositionsBlock(CellBlock &block, double delta_t) {
    block.forEachCell([delta_t](Cell &cell) {
        updatePositionsCell(cell, delta_t);
    });
}

void updateVelocitiesBlock(CellBlock &block, double delta_t) {
    block.forEachCell([delta_t](Cell &cell) {
        updateVelocitiesCell(cell, delta_t);
    });
}

