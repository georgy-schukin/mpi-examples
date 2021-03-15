#pragma once

#include <vector>

class Slicer {
public:
    static std::vector<int> getSlices(int size, int num);
    static std::vector<int> getShifts(int size, int num);
    static std::vector<int> getShifts(const std::vector<int> &sizes);
    static int shiftIndexFor(const std::vector<int> &shifts, int value);
};
