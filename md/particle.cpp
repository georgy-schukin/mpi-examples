#include "particle.h"

bool Particle::isIn(const Extent3 &extent) const {
    return extent.isIn(0, pos.x) &&
           extent.isIn(1, pos.y) &&
           extent.isIn(2, pos.z);
}

