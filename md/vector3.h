#pragma once

#include "comm/serialize/serializable.h"

class Vector3 {
public:
    Vector3() {}
    Vector3(double x, double y, double z) :
        x(x), y(y), z(z) {
    }

    inline void null() {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }

    inline void add(const Vector3 &v) {
        x += v.x;
        y += v.y;
        z += v.z;
    }

    inline void add(double xx, double yy, double zz) {
        x += xx;
        y += yy;
        z += zz;
    }

    inline void sub(const Vector3 &v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
    }

    inline void sub(double xx, double yy, double zz) {
        x -= xx;
        y -= yy;
        z -= zz;
    }

public:
    double x = 0.0, y = 0.0, z = 0.0;
};

ddl::Writable& operator<< (ddl::Writable &w, const Vector3 &v);
ddl::Readable& operator>> (ddl::Readable &r, Vector3 &v);
