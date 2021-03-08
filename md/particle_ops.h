#pragma once

#include "particle.h"
#include "extent.h"
#include "index.h"

#include <vector>

struct MDParams;
struct MeshParams;

void randomPosition(std::vector<Particle> &particles, const Extent3 &extent);
void randomPositionEllipse(std::vector<Particle> &particles, const Extent3 &extent);
void regularMeshPosition(std::vector<Particle> &particles, const Extent3 &extent, const Vector3 &step);
Index3 regularMeshSize(const Extent3 &extent, const Vector3 &step);

void uniformMass(std::vector<Particle> &particles, double mass);
void randomMass(std::vector<Particle> &particles, double min, double max);

void uniformVelocity(std::vector<Particle> &particles, const Vector3 &velocity);
void randomVelocity(std::vector<Particle> &particles, double scale = 1.0);

std::vector<Particle> initRandom(int num_of_particles, const MDParams &md_params, const MeshParams &mesh_params);
std::vector<Particle> initRandomEllipse(int num_of_particles, const MDParams &md_params, const MeshParams &mesh_params);
std::vector<Particle> initCollision(const MDParams &md_params, const MeshParams &mesh_params);
std::vector<Particle> initExplosion(int num_of_particles, const MDParams &md_params, const MeshParams &mesh_params);

std::vector<Particle> initParticles(int num_of_particles, const MDParams &md_params, const MeshParams &mesh_params);
