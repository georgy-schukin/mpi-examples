#include "slicer.h"

std::vector<int> Slicer::getSlices(int size, int num) {
    const int block_size = size / num;
    const int remaining = size % num;
    std::vector<int> slice_sizes;
    for (int i = 0; i < num; i++) {
        const int slice_size = (i < remaining ? block_size + 1 : block_size);
        slice_sizes.push_back(slice_size);
    }
    return slice_sizes;
}

std::vector<int> Slicer::getShifts(int size, int num) {
    return getShifts(getSlices(size, num));
}

std::vector<int> Slicer::getShifts(const std::vector<int> &sizes) {
    std::vector<int> shifts = {0};
    for (auto slice_size: sizes) {
        const auto shift = shifts.back() + slice_size;
        shifts.push_back(shift);
    }
    return shifts;
}

int Slicer::shiftIndexFor(const std::vector<int> &shifts, int value) {
    for (int i = 0; i < (int)(shifts.size() - 1); i++)   {
        if (value >= shifts[i] && value < shifts[i + 1]) {
            return i;
        }
    }
    return -1;
}
