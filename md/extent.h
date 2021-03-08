#pragma once

#include "vector3.h"

#include <array>
#include <cstddef>

template <typename T, size_t Dims>
class Extent {
public:
    Extent() {}

    const T& start(int dim) const {
        return _start[dim];
    }

    T& start(int dim) {
        return _start[dim];
    }

    const T& end(int dim) const {
        return _end[dim];
    }

    T& end(int dim) {
        return _end[dim];
    }

    T size(int dim) const {
        return end(dim) - start(dim);
    }

    bool isIn(int dim, const T &value) const {
        return value >= start(dim) && value < end(dim);
    }

    bool isOut(int dim, const T &value) const {
        return value < start(dim) || value >= end(dim);
    }

    Extent scaledFromCenter(double scale) const {
        Extent scaled;
        for (int i = 0; i < Dims; i++) {
            const auto center = start(i) + size(i) * 0.5;
            scaled.start(i) = center - size(i) * 0.5 * scale;
            scaled.end(i) = center + size(i) * 0.5 * scale;
        }
        return scaled;
    }

    Extent scaledFromCenter(const Vector3 &scale) const {
        Extent scaled;
        std::array<double, 3> sc = {scale.x, scale.y, scale.z};
        for (int i = 0; i < Dims; i++) {
            const auto center = start(i) + size(i) * 0.5;
            scaled.start(i) = center - size(i) * 0.5 * sc[i];
            scaled.end(i) = center + size(i) * 0.5 * sc[i];
        }
        return scaled;
    }

    Extent shifted(const Vector3 &shift) {
        std::array<double, 3> sf = {shift.x, shift.y, shift.z};
        Extent shifted_extent;
        for (int i = 0; i < 3; i++) {
            shifted_extent.start(i) = start(i) + sf[i];
            shifted_extent.end(i) = end(i) + sf[i];
        }
        return shifted_extent;
    }

    Vector3 getCenter() const {
        return Vector3(start(0) + size(0) * 0.5,
                       start(1) + size(1) * 0.5,
                       start(2) + size(2) * 0.5);
    }

private:
    std::array<T, Dims> _start;
    std::array<T, Dims> _end;
};
