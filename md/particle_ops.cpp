#include "particle_ops.h"
#include "md.h"

#include <random>
#include <cmath>

void randomPosition(std::vector<Particle> &particles, const Extent3 &extent) {
    std::mt19937 gen;
    std::uniform_real_distribution<double> distr(0.0, 1.0);
    for (auto &particle: particles) {
        const auto px = extent.start(0) + distr(gen) * extent.size(0);
        const auto py = extent.start(1) + distr(gen) * extent.size(1);
        const auto pz = extent.start(2) + distr(gen) * extent.size(2);
        particle.pos = Vector3(px, py, pz);
    }
}

void randomPositionEllipse(std::vector<Particle> &particles, const Extent3 &extent) {
    std::mt19937 gen;
    std::uniform_real_distribution<double> distr(0.0, 1.0);
    const double a = extent.size(0) * 0.5;
    const double b = extent.size(1) * 0.5;
    const double c = extent.size(2) * 0.5;
    for (auto &particle: particles) {
        double px = 0, py = 0, pz = 0;
        double x = 0, y = 0, z = 0;
        do {
            px = extent.start(0) + distr(gen) * extent.size(0);
            py = extent.start(1) + distr(gen) * extent.size(1);
            pz = extent.start(2) + distr(gen) * extent.size(2);
            x = px - (extent.start(0) + a);
            y = py - (extent.start(1) + b);
            z = pz - (extent.start(2) + c);
        } while (x * x / (a * a) + y * y / (b * b) + z * z / (c * c) > 1.0);
        particle.pos = Vector3(px, py, pz);
    }
}

void regularMeshPosition(std::vector<Particle> &particles, const Extent3 &extent, const Vector3 &step) {
    auto size = regularMeshSize(extent, step);
    size_t index = 0;
    for (int i = 0; i < size[0]; i++)
    for (int j = 0; j < size[1]; j++)
    for (int k = 0; k < size[2]; k++) {
        if (index >= particles.size()) {
            return;
        }
        auto &particle = particles[index];
        particle.pos = Vector3(extent.start(0) + i * step.x,
                               extent.start(1) + j * step.y,
                               extent.start(2) + k * step.z);
        index++;
    }
}

Index3 regularMeshSize(const Extent3 &extent, const Vector3 &step) {
    return Index3 {static_cast<int>(extent.size(0) / step.x) + 1,
                   static_cast<int>(extent.size(1) / step.y) + 1,
                   static_cast<int>(extent.size(2) / step.z) + 1};
}

void uniformMass(std::vector<Particle> &particles, double mass) {
    for (auto &particle: particles) {
        particle.mass = mass;
    }
}

void randomMass(std::vector<Particle> &particles, double min, double max) {
    std::mt19937 gen;
    std::uniform_real_distribution<double> distr(min, std::nextafter(max, std::numeric_limits<double>::max()));
    for (auto &particle: particles) {
        particle.mass = distr(gen);
    }
}

void uniformVelocity(std::vector<Particle> &particles, const Vector3 &velocity) {
    for (auto &particle: particles) {
        particle.velocity = velocity;
    }
}

void randomVelocity(std::vector<Particle> &particles, double scale) {
    std::mt19937 gen;
    std::uniform_real_distribution<double> distr(-1.0, std::nextafter(1.0, std::numeric_limits<double>::max()));
    for (auto &particle: particles) {
        const auto vx = distr(gen);
        const auto vy = distr(gen);
        const auto vz = distr(gen);
        const auto len = std::sqrt(vx * vx + vy * vy + vz * vz) + 1e-12;
        const auto s = scale / len;
        particle.velocity = Vector3(vx * s, vy * s, vz * s);
    }
}

std::vector<Particle> initRandom(int num_of_particles, const MDParams &md_params, const MeshParams &mesh_params) {
    std::vector<Particle> particles(static_cast<size_t>(num_of_particles));
    randomPosition(particles, mesh_params.getExtent());
    randomMass(particles, md_params.min_mass, md_params.max_mass);
    randomVelocity(particles, md_params.init_speed);
    return particles;
}

std::vector<Particle> initRandomEllipse(int num_of_particles, const MDParams &md_params, const MeshParams &mesh_params) {
    std::vector<Particle> particles(static_cast<size_t>(num_of_particles));
    randomPositionEllipse(particles, mesh_params.getExtent());
    randomMass(particles, md_params.min_mass, md_params.max_mass);
    randomVelocity(particles, md_params.init_speed);
    return particles;
}

std::vector<Particle> initCollision(const MDParams &md_params, const MeshParams &mesh_params) {
    const auto step = std::pow(2.0, 1.0/6.0) * md_params.sigma;
    Vector3 steps {step, step, step};
    const auto extent = mesh_params.getExtent();
    const auto big_body_extent = extent.
            scaledFromCenter(Vector3 {0.8, 0.8, 0.3}).
            shifted(Vector3 {0, 0, -0.125 * extent.size(2)});
    const auto small_body_extent = mesh_params.getExtent().
            scaledFromCenter(Vector3 {0.4, 0.4, 0.25}).
            shifted(Vector3 {0, 0, 0.3 * extent.size(2)});
    const auto big_body_size = regularMeshSize(big_body_extent, steps);
    const auto small_body_size = regularMeshSize(small_body_extent, steps);
    std::vector<Particle> big_body(static_cast<size_t>(big_body_size[0] * big_body_size[1] * big_body_size[2]));
    std::vector<Particle> small_body(static_cast<size_t>(small_body_size[0] * small_body_size[1] * small_body_size[2]));
    regularMeshPosition(big_body, big_body_extent, steps);
    regularMeshPosition(small_body, small_body_extent, steps);
    uniformMass(big_body, md_params.max_mass);
    uniformMass(small_body, md_params.min_mass);
    uniformVelocity(big_body, Vector3 {});
    uniformVelocity(small_body, Vector3 {0.0, 0.0, -md_params.init_speed});
    std::vector<Particle> particles;
    particles.insert(particles.end(), big_body.begin(), big_body.end());
    particles.insert(particles.end(), small_body.begin(), small_body.end());
    return particles;
}

std::vector<Particle> initExplosion(int num_of_particles, const MDParams &md_params, const MeshParams &mesh_params) {
    std::vector<Particle> particles(static_cast<size_t>(num_of_particles));
    const auto extent = mesh_params.getExtent().scaledFromCenter(0.6);
    randomPositionEllipse(particles, extent);
    const auto center = extent.getCenter();
    const auto max_distance = std::max(extent.size(0), std::max(extent.size(1),extent.size(2))) * 0.5;
    for (auto &particle: particles) {
        const auto &p = particle.pos;
        const Vector3 v {p.x - center.x, p.y - center.y, p.z - center.z};
        const auto distance_squared = v.x * v.x + v.y * v.y + v.z * v.z + 1e-12;
        const auto t = distance_squared / (max_distance * max_distance);
        particle.mass = t * md_params.min_mass + (1.0 - t) * md_params.max_mass;
        const auto vel_c = md_params.init_speed * (1.0 - t) / std::sqrt(distance_squared);
        particle.velocity = Vector3(vel_c * v.x, vel_c * v.y, vel_c * v.z);
    }
    return particles;
}

std::vector<Particle> initParticles(int num_of_particles, const MDParams &md_params, const MeshParams &mesh_params) {
    switch (md_params.particles_init) {
        case MDParams::INIT_RANDOM: return initRandom(num_of_particles, md_params, mesh_params);
        case MDParams::INIT_RANDOM_ELLIPSE: return initRandomEllipse(num_of_particles, md_params, mesh_params);
        case MDParams::INIT_COLLISION: return initCollision(md_params, mesh_params);
        case MDParams::INIT_EXPLOSION: return initExplosion(num_of_particles, md_params, mesh_params);
        default: return initRandom(num_of_particles, md_params, mesh_params);
    }
}


