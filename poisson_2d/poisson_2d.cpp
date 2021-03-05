#include "poisson_2d.h"

namespace psn2d {

void initData(std::vector<double> &data, const DataParams &params, const PoissonParams &p_params, Func3D init) {
    #define EL(data,i,j,k) data[((i)*params.total_y*params.total_z + (j)*params.total_z + (k))]
	
    for (int i = params.shift_x; i < params.dim_x + params.shift_x; i++)
    for (int j = params.shift_y; j < params.dim_y + params.shift_y; j++)
    for (int k = 0; k < p_params.grid_size_z; k++) {
    	if ((i == 0) || (i == p_params.grid_size_x - 1) || 
    	    (j == 0) || (j == p_params.grid_size_y - 1) || 
    	    (k == 0) || (k == p_params.grid_size_z - 1)) {
    	    EL(data, i - params.shift_x + params.shadow_x[0], 
    	             j - params.shift_y + params.shadow_y[0], 
    	             k) = init(i * p_params.step_x, j * p_params.step_y, k * p_params.step_z); // bound values
    	}
    }
    	
    #undef EL
}

void initRight(std::vector<double> &right, const DataParams &params, const PoissonParams &p_params, Func3D func) {
    #define EL(data,i,j,k) data[((i)*params.total_y*params.total_z + (j)*params.total_z + (k))]
		
    for (int i = params.shift_x; i < params.dim_x + params.shift_x; i++)
    for (int j = params.shift_y; j < params.dim_y + params.shift_y; j++)
    for (int k = 0; k < p_params.grid_size_z; k++) {
        EL(right, i - params.shift_x + params.shadow_x[0], 
    	          j - params.shift_y + params.shadow_y[0], 
    	          k) = func(i * p_params.step_x, j * p_params.step_y, k * p_params.step_z);
    }
    	
    #undef EL
}

void computeData(std::vector<double> &data, const std::vector<double> &prev_data, const std::vector<double> &right, const DataParams &params, const PoissonParams &p_params) {
    #define EL(data,i,j,k) data[((i)*params.total_y*params.total_z + (j)*params.total_z + (k))]
    const double OWX = p_params.step_x * p_params.step_x;
    const double OWY = p_params.step_y * p_params.step_y;
    const double OWZ = p_params.step_z * p_params.step_z;
    const double OWX1 = 1.0/OWX;
    const double OWY1 = 1.0/OWY;
    const double OWZ1 = 1.0/OWZ;
    const double C = 2.0*OWX1 + 2.0*OWY1 + 2.0*OWZ1;
    const double C1 = 1.0/C;
    for (int i = 1; i < params.total_x - 1; i++)
    for (int j = 1; j < params.total_y - 1; j++)
	for (int k = 1; k < params.total_z - 1; k++) {        	    	    	    	    
	    const double Fi = (EL(prev_data,i - 1,j,k) + EL(prev_data,i + 1,j,k)) * OWX1;
	    const double Fj = (EL(prev_data,i,j - 1,k) + EL(prev_data,i,j + 1,k)) * OWY1;
	    const double Fk = (EL(prev_data,i,j,k - 1) + EL(prev_data,i,j,k + 1)) * OWZ1;
	    	     	    	     
	    const double R = EL(right,i,j,k);
	    	     
	    EL(data,i,j,k) = (Fi + Fj + Fk - R) * C1;
    }
    #undef EL    	
}

}
