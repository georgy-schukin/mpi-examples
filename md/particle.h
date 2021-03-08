#pragma once

#include "vector3.h"
#include "extent.h"

using Extent3 = Extent<double, 3>;

class Particle {
public:
    Particle() {}

    bool isIn(const Extent3 &extent) const;

public:
    Vector3 pos;
    Vector3 velocity;
    Vector3 force;
    Vector3 force_prev;
    double mass = 1.0;
};

