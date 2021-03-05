#pragma once

#include <vector>
#include <functional>

namespace psn2d {

const double X = 1.0; // size of data area by X
const double Y = 1.0; // size of data area by Y
const double Z = 1.0; // size of data area by Z

const int REPEAT = 200; // repeat computations this times

const double A = 1.0;
const double EPS = 1e-6;

struct PoissonParams {
    int grid_size_x;
    int grid_size_y;
    int grid_size_z;
    double step_x;
    double step_y;
    double step_z;
};

// Parameters of data fragment (block)
struct DataParams {
    int rank_x;
    int rank_y;
    int dim_x;
    int dim_y;
    int shift_x;
    int shift_y;
    int total_x;
    int total_y;
    int total_z;
    int shadow_x[2];
    int shadow_y[2];
};

using Func3D = std::function<double(double, double, double)>;

void initData(std::vector<double> &data, const DataParams &params, const PoissonParams &p_params, Func3D init);
void initRight(std::vector<double> &right, const DataParams &params, const PoissonParams &p_params, Func3D func);
void computeData(std::vector<double> &data, const std::vector<double> &prev_data, const std::vector<double> &right, const DataParams &params, const PoissonParams &p_params);

}