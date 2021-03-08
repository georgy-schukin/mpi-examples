#include "cell.h"
#include "md.h"

Extent3 Cell::getExtent(const MeshParams &m_params) const {
    Extent3 extent;
    for (int i = 0; i < 3; i++) {
        extent.start(i) = m_params.origin[i] + index[i] * m_params.step[i];
        extent.end(i) = extent.start(i) + m_params.step[i];
    }
    return extent;
}

std::vector<Particle> Cell::extractOutParticles(const MeshParams &m_params) {
    const auto extent = getExtent(m_params);    
    std::vector<Particle> in_particles, out_particles;
    std::lock_guard<std::mutex> lock(p_mutex);
    for (const auto &particle: particles) {
        if (particle.isIn(extent)) {
            in_particles.push_back(particle);
        } else {
            out_particles.push_back(particle);
        }
    }        
    particles = std::move(in_particles);
    return out_particles;
}

void Cell::addParticles(const std::vector<Particle> &new_particles) {
    std::lock_guard<std::mutex> lock(p_mutex);
    particles.insert(particles.end(), new_particles.begin(), new_particles.end());
}

void Cell::clearParticles() {
    std::lock_guard<std::mutex> lock(p_mutex);
    particles.clear();
}

