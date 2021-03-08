#pragma once

#include "particle.h"
#include "extent.h"
#include "index.h"

#include <vector>
#include <mutex>
#include <ostream>

struct MeshParams;

class Cell {
public:
    Cell() {}

    Cell(const Cell &c) :
        index(c.index), particles(c.particles), neigh_indices_type(c.neigh_indices_type) {
    }

    Cell(Cell &&c) :
        index(c.index), particles(std::move(c.particles)), neigh_indices_type(c.neigh_indices_type) {
    }

    Cell& operator=(const Cell &c) {        
        index = c.index;
        particles = c.particles;
        neigh_indices_type = c.neigh_indices_type;
        return *this;
    }

    Cell& operator=(Cell &&c) {        
        index = c.index;
        particles = std::move(c.particles);
        neigh_indices_type = c.neigh_indices_type;
        return *this;
    }

    std::vector<Particle>& getParticles() {
        return particles;
    }

    const std::vector<Particle>& getParticles() const {
        return particles;
    }

    const Index3& getIndex() const {
        return index;
    }

    void setIndex(const Index3 &index) {
        this->index = index;
    }

    void setNeighborIndicesType(int type) {
        neigh_indices_type = type;
    }

    int getNeighborIndicesType() const {
        return neigh_indices_type;
    }

    Extent3 getExtent(const MeshParams &m_params) const;

    std::vector<Particle> extractOutParticles(const MeshParams &m_params);
    void addParticles(const std::vector<Particle> &new_particles);
    void clearParticles();

    template <typename Proc, typename... Args>
    void forEachParticle(Proc proc, Args... args) {
        for (auto &particle: getParticles()) {
            proc(particle, args...);
        }
    }

    template <typename Proc, typename... Args>
    void forEachParticle(Proc proc, Args... args) const {
        for (const auto &particle: getParticles()) {
            proc(particle, args...);
        }
    }

private:
    Index3 index; // cell's index in 3D mesh
    std::vector<Particle> particles;
    int neigh_indices_type = 0;
    std::mutex p_mutex;
};
