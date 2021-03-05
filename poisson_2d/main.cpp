/**
** 3D Puasson with 2D decomposition by nodes.
*/

#include <mpi.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <sstream>

#include "poisson_2d.h"

using namespace psn2d;

double Fresh(const double &x, const double &y, const double &z) {
    return x + y + z;
}

double Ro(const double &x, const double &y, const double &z) {
    return 0; //-psn2d::A*(x + y + z);
}

void slice(const int size, const int num, std::vector<int> &sizes, std::vector<int> &shifts) {
    for (int i = 0; i < num; i++)
	sizes[i] = size/num + (i < size % num ? 1 : 0);
	
    shifts[0] = 0;
    for (int i = 1; i < num; i++)
	shifts[i] = shifts[i - 1] + sizes[i - 1];	
}

int fromBoolean(const bool value) {
    return value ? 1 : 0;
}

int main(int argc, char** argv) {		
	MPI_Comm cart_comm;
	int rank, size, size_x, size_y;
	
	MPI_Init(&argc, &argv);
	MPI_Request req[8];
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	int grid_size = (argc > 1 ? atoi(argv[1]) : 100);
	int max_iter = (argc > 2 ? atoi(argv[2]) : 10);
	bool scale = (argc > 3 ? atoi(argv[3]) : 0);
	
	// Find optimal sizes of process grid
	for (int i = int(sqrtf((float)size)); i >= 1; i--) {
	    if (size % i == 0) {
		  size_x = i;
		  size_y = size / i;
		  break;
	    }  
	}		
		    		
	int dims[2] = {size_x, size_y};
	int periods[2] = {0, 0};
	
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);
	
	MPI_Comm_rank(cart_comm, &rank);
	
	int coords[2];
	MPI_Cart_coords(cart_comm, rank, 2, coords); // get process indices by X and Y
	
	const int rank_x = coords[0];
	const int rank_y = coords[1]; 
			
	PoissonParams p_params;
	p_params.grid_size_x = scale ? grid_size * size_x : grid_size;
	p_params.grid_size_y = scale ? grid_size * size_y : grid_size;
	p_params.grid_size_z = grid_size;	
	p_params.step_x = X / p_params.grid_size_x;
	p_params.step_y = X / p_params.grid_size_y;	
	p_params.step_z = X / p_params.grid_size_z;	
	
	std::vector<int> dft_x(size_x), sft_x(size_x);
	std::vector<int> dft_y(size_y), sft_y(size_y);	
	
	slice(p_params.grid_size_x, size_x, dft_x, sft_x);
	slice(p_params.grid_size_y, size_y, dft_y, sft_y);
	
	DataParams params;
	params.rank_x = rank_x;
	params.rank_y = rank_y;
	params.dim_x = dft_x[rank_x]; // num of elements by X in node data piece
	params.dim_y = dft_y[rank_y];
	params.shift_x = sft_x[rank_x]; // shift of elements indices from 0 by X in node piece
	params.shift_y = sft_y[rank_y];
	params.shadow_x[0] = fromBoolean(rank_x > 0);
	params.shadow_x[1] = fromBoolean(rank_x < size_x - 1);
	params.shadow_y[0] = fromBoolean(rank_y > 0);
	params.shadow_y[1] = fromBoolean(rank_y < size_y - 1);	
	params.total_x = params.dim_x + params.shadow_x[0] + params.shadow_x[1]; // total num of elements by X (shadows included) in node piece
	params.total_y = params.dim_y + params.shadow_y[0] + params.shadow_y[1];
	params.total_z = p_params.grid_size_z;
	
	//printf("%d: %d %d, data size %d %d, total %d %d\n", rank, rank_x, rank_y, params.dim_x, params.dim_y, params.total_x, params.total_y);		
		
	std::vector<double> data(params.total_x * params.total_y * params.total_z); // Data
	std::vector<double> data_new(params.total_x * params.total_y * params.total_z); // New data
	std::vector<double> right(params.total_x * params.total_y * params.total_z); // Right part
	
	MPI_Datatype x_type, y_type; // type for shadow exchange by X and Y
	MPI_Type_vector(params.total_y, params.total_z, 0, MPI_DOUBLE, &x_type); 
	MPI_Type_vector(params.total_x, params.total_z, params.total_y * params.total_z, MPI_DOUBLE, &y_type);
	MPI_Type_commit(&x_type);
	MPI_Type_commit(&y_type);							    
  	 
	const double start_time = MPI_Wtime();
	
	initData(data, params, p_params, Fresh);
	initRight(right, params, p_params, Ro);
	
	//double N = 0.0, MAK = 0.0, Fa = 0.0;
	
	int it = 0;
	double global_max = 0.0;
		
	int neigh_x[2], neigh_y[2];
	MPI_Cart_shift(cart_comm, 0, 1, &neigh_x[0], &neigh_x[1]); // by X: left nd right neighbour
	MPI_Cart_shift(cart_comm, 1, 1, &neigh_y[0], &neigh_y[1]); // by Y: top and bottom neighbour
	
	//printf("%d: x %d %d, y %d %d\n", rank, neigh_x[0], neigh_x[1], neigh_y[0], neigh_y[1]);       
	#define EL(data,i,j,k) data[(i)*params.total_y*params.total_z + (j)*params.total_z + (k)]
	
	double comp_time = 0, update_time = 0, copy_time = 0;
	
      do { 	    	    
        double t_us = MPI_Wtime();
	    if (params.shadow_x[1]) // update right shadow by x
	    {
        	    MPI_Isend(&EL(data,params.total_x - 2,0,0), 1, x_type, neigh_x[1], 1, cart_comm, &req[0]); // send data to shadow
        	    MPI_Irecv(&EL(data,params.total_x - 1,0,0), 1, x_type, neigh_x[1], 2, cart_comm, &req[1]); // recv data in shadow
    	    }
    	    
    	    if (params.shadow_x[0]) // update left shadow by x
    	    {
    	        MPI_Isend(&EL(data,1,0,0), 1, x_type, neigh_x[0], 2, cart_comm, &req[2]); // send to shadow
        	    MPI_Irecv(&EL(data,0,0,0), 1, x_type, neigh_x[0], 1, cart_comm, &req[3]); // recv in shadow
    	    }

	    if (params.shadow_y[1]) // update top shadow by y
    	    {
        	    MPI_Isend(&EL(data,0,params.total_y - 2,0), 1, y_type, neigh_y[1], 3, cart_comm, &req[4]); // send data to shadow
        	    MPI_Irecv(&EL(data,0,params.total_y - 1,0), 1, y_type, neigh_y[1], 4, cart_comm, &req[5]); // recv data in shadow
    	    }
    	    
    	    if (params.shadow_y[0]) // update bottom shadow by y
    	    {
    	        MPI_Isend(&EL(data,0,1,0), 1, y_type, neigh_y[0], 4, cart_comm, &req[6]); // send to shadow
        	    MPI_Irecv(&EL(data,0,0,0), 1, y_type, neigh_y[0], 3, cart_comm, &req[7]); // recv in shadow
    	    }
    	    
    	    if (params.shadow_x[1]) {
    		    MPI_Wait(&req[0], MPI_STATUS_IGNORE);
    		    MPI_Wait(&req[1], MPI_STATUS_IGNORE);
    	    }
    	    
    	    if (params.shadow_x[0]) {
    		    MPI_Wait(&req[2], MPI_STATUS_IGNORE);
    		    MPI_Wait(&req[3], MPI_STATUS_IGNORE);
    	    }

    	    if (params.shadow_y[1]) {
    		    MPI_Wait(&req[4], MPI_STATUS_IGNORE);
    		    MPI_Wait(&req[5], MPI_STATUS_IGNORE);
    	    }
    	    
    	    if (params.shadow_y[0]) {
    		    MPI_Wait(&req[6], MPI_STATUS_IGNORE);
    		    MPI_Wait(&req[7], MPI_STATUS_IGNORE);
    	    }
    	double t_ue = MPI_Wtime();  
    	update_time += (t_ue - t_us);
    	    
	    double max = 0.0;	  
	    
	    double t_cs = MPI_Wtime();
	    computeData(data_new, data, right, params, p_params);
	    double t_ce = MPI_Wtime();  
    	comp_time += (t_ce - t_cs);
        	
       	//global_max = 0.0;
    	//MPI_Allreduce(&max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);	    	    
    	    
        double t_cps = MPI_Wtime();
    	data = data_new;  
    	double t_cpe = MPI_Wtime();  
    	copy_time += (t_cpe - t_cps);
        	
    	it++;              
	}
	while(/*(global_max > EPS) &&*/ (it < max_iter)); 
	
        const double end_time = MPI_Wtime();
        
        double local_time = end_time - start_time;
        double max_time = 0;
        
        MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);	                                  
        
        if (rank == 0) {
        	printf("Nodes: %d\n", size);
	        printf("Grid %d %d\n", size_x, size_y);
    	    printf("Size: %d %d %d\n", p_params.grid_size_x, p_params.grid_size_y, p_params.grid_size_z);
    	    printf("Iterations: %d\n", it);        
    	    printf("TIME: %.5f\n", max_time);
        }
        
      printf("Node %d: comp: %.5f, update: %.5f, copy: %.5f\n", rank, comp_time, update_time, copy_time);        
                  
	// Нахождение максимального расхождения полученного приближенного решения
	// и точного решения                  
/*	double max = 0.0;                 
	for (int i = 1; i < params.total_x - 1; i++)
    	for (int j = 1; j < params.total_y - 1; j++)
    	for (int k = 1; k < SZ - 1; k++) {
            double F1 = fabs(F(L1,i,j,k) - Fresh((i - params.shadow_x[0] + params.shift_x)*HX, (j - params.shadow_y[0] + params.shift_y)*HY, k*HZ));
            if(F1 > max) max = F1;                    
    	}
	
	global_max = 0;
	MPI_Allreduce(&max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);	                                  
    
	if(rank == 0)
	    printf("Max differ = %e\n", global_max);*/
        
        MPI_Type_free(&x_type);
        MPI_Type_free(&y_type);
	MPI_Finalize();
	
	return 0;          
} 
 