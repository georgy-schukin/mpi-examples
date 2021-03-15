#include <mpi.h>

#include "mesh3d.h"
#include "md.h"
#include "output.h"
#include "neighbor.h"
#include "particle_ops.h"

#include <iostream>
#include <string>
#include <sstream>
#include <chrono>
#include <stdexcept>
#include <fstream>
#include <cmath>

int intArg(int arg_num, int default_value, int argc, char **argv) {
    return (argc > arg_num ? std::stoi(argv[arg_num]) : default_value);
}

int toRank(const Index3 &index, const MeshDistribution &distr) {
    return index[0] * static_cast<int>(distr.numOfNodes(1)) * static_cast<int>(distr.numOfNodes(2)) +
           index[1] * static_cast<int>(distr.numOfNodes(2)) +
           index[2];
}

double inExtent(double value, double start, double end) {
    if (value < start) {
        return start + std::fmod(start - value, end - start);
    } else if (value >= end) {
        return start + std::fmod(value - start, end - start);
    } else {
        return value;
    }
}

Index3 getBlockIndex(const Vector3 &pos, const MeshParams &mesh_params, const MeshDecomosition &mesh_decomp) {
    Index3 index;
    std::array<double, 3> p {pos.x, pos.y, pos.z};
    for (size_t i = 0; i < 3; i++) {
        const auto cell_index = int((p[i] - mesh_params.origin[i]) / mesh_params.step[i]);
        index[i] = mesh_decomp.blockIndex(static_cast<int>(i), cell_index);
    }
    return index;
}

void distributeParticles(Mesh3D &mesh, const std::vector<Particle> &particles, const MeshParams &mesh_params,
                           const MeshDecomosition &decomp, const std::map<Index3, int> &block_distr) {
    static auto particle_add = [mesh_params](CellBlock &block, const std::vector<Particle> &particles) {
        block.addParticles(particles, mesh_params);
    };
    static auto add_particles = mesh.makeGetter<decltype(particle_add), void, const std::vector<Particle>&>(particle_add);

    std::map<Index3, std::vector<Particle>> particles_for_blocks;

    for (const auto &particle: particles) {
        particles_for_blocks[getBlockIndex(particle.pos, mesh_params, decomp)].push_back(particle);
    }

    std::vector<std::future<void>> futures;
    for (auto &p: particles_for_blocks) {
        auto it = block_distr.find(p.first);
        if (it != block_distr.end()) {
            futures.push_back(add_particles(it->second, it->first, p.second));
        }
    }
    // Wait when remote operations are completed.
    for (auto &f: futures) {
        f.get();
    }
}

template <typename Proc>
void initParticles(Proc proc, Mesh3D &mesh, int my_rank, const MeshParams &mesh_params,
                   const MeshDecomosition &decomp, const std::map<Index3, int> &block_distr) {
    std::vector<Particle> particles;
    if (my_rank == 0) {
        particles = proc();
    }
    distributeParticles(mesh, particles, mesh_params, decomp, block_distr);
}

void moveParticles(Mesh3D &mesh, const MeshParams &mesh_params,
                   const MeshDecomosition &decomp, const std::map<Index3, int> &block_distr) {
    const auto extent = mesh_params.getExtent();
    const bool account_periodic = mesh_params.periodic[0] || mesh_params.periodic[1] || mesh_params.periodic[2];

    std::vector<Particle> moved_particles;
    for (auto &p: mesh.localContents()) {
        auto out_particles = p.second->extractAndMoveOutParticles(mesh_params);
        if (account_periodic) {
            // Change particle positions to account periodicity.
            for (auto &particle: out_particles) {
                if (!particle.isIn(extent)) {
                    particle.pos.x = mesh_params.periodic[0] ? inExtent(particle.pos.x, extent.start(0), extent.end(0)) : particle.pos.x;
                    particle.pos.y = mesh_params.periodic[1] ? inExtent(particle.pos.y, extent.start(1), extent.end(1)) : particle.pos.y;
                    particle.pos.z = mesh_params.periodic[2] ? inExtent(particle.pos.z, extent.start(2), extent.end(2)) : particle.pos.z;
                }
            }
        }
        moved_particles.insert(moved_particles.end(), out_particles.begin(), out_particles.end());
    }

    distributeParticles(mesh, moved_particles, mesh_params, decomp, block_distr);
}

enum Dimension {
    DIM_X = 0,
    DIM_Y,
    DIM_Z
};

struct TransferInput {
    int dimension;
    int src1, src2, src3;
    int dst1, dst2, dst3;
    int size1, size2;

    TransferInput() {}
    TransferInput(int dim, int s1, int s2, int s3,
                  int d1, int d2, int d3, int sz1, int sz2) :
        dimension(dim), src1(s1), src2(s2), src3(s3),
        dst1(d1), dst2(d2), dst3(d3), size1(sz1), size2(sz2) {
    }
};

std::map<Index3, TransferInput> getTransferInputs(const CellBlockParams &params, const NeighborIndices &neigh_indices) {
    const auto &shadow = params.shadow_start;
    const auto size = params.size;
    static const int NEXT = 1;
    static const int PREV = -1;
    static const int START = 0;
    static const int END = -1;
    static const int S_PREV = END;
    static const int S_NEXT = START;
    static const int D_PREV = START;
    static const int D_NEXT = END;
    std::map<Index3, TransferInput> all_transfers = {
        {{PREV, 0, 0}, TransferInput(DIM_X, S_PREV, 0, 0, D_PREV, shadow[1], shadow[2], size[1], size[2])}, // front X
        {{NEXT, 0, 0}, TransferInput(DIM_X, S_NEXT, 0, 0, D_NEXT, shadow[1], shadow[2], size[1], size[2])}, // back X
        {{0, PREV, 0}, TransferInput(DIM_Y, 0, S_PREV, 0, shadow[0], D_PREV, shadow[2], size[0], size[2])}, // front Y
        {{0, NEXT, 0}, TransferInput(DIM_Y, 0, S_NEXT, 0, shadow[0], D_NEXT, shadow[2], size[0], size[2])}, // back Y
        {{0, 0, PREV}, TransferInput(DIM_Z, 0, 0, S_PREV, shadow[0], shadow[1], D_PREV, size[0], size[1])}, // front Z
        {{0, 0, NEXT}, TransferInput(DIM_Z, 0, 0, S_NEXT, shadow[0], shadow[1], D_NEXT, size[0], size[1])}, // back Z
        {{PREV, PREV, 0}, TransferInput(DIM_X, S_PREV, S_PREV, 0, D_PREV, D_PREV, shadow[2], 1, size[2])}, // sides XY (Z changes)
        {{PREV, NEXT, 0}, TransferInput(DIM_X, S_PREV, S_NEXT, 0, D_PREV, D_NEXT, shadow[2], 1, size[2])},
        {{NEXT, PREV, 0}, TransferInput(DIM_X, S_NEXT, S_PREV, 0, D_NEXT, D_PREV, shadow[2], 1, size[2])},
        {{NEXT, NEXT, 0}, TransferInput(DIM_X, S_NEXT, S_NEXT, 0, D_NEXT, D_NEXT, shadow[2], 1, size[2])},
        {{PREV, 0, PREV}, TransferInput(DIM_X, S_PREV, 0, S_PREV, D_PREV, shadow[1], D_PREV, size[1], 1)}, // sides XZ (Y changes)
        {{PREV, 0, NEXT}, TransferInput(DIM_X, S_PREV, 0, S_NEXT, D_PREV, shadow[1], D_NEXT, size[1], 1)},
        {{NEXT, 0, PREV}, TransferInput(DIM_X, S_NEXT, 0, S_PREV, D_NEXT, shadow[1], D_PREV, size[1], 1)},
        {{NEXT, 0, NEXT}, TransferInput(DIM_X, S_NEXT, 0, S_NEXT, D_NEXT, shadow[1], D_NEXT, size[1], 1)},
        {{0, PREV, PREV}, TransferInput(DIM_Y, 0, S_PREV, S_PREV, shadow[0], D_PREV, D_PREV, size[0], 1)}, // sides YZ (X changes)
        {{0, PREV, NEXT}, TransferInput(DIM_Y, 0, S_PREV, S_NEXT, shadow[0], D_PREV, D_NEXT, size[0], 1)},
        {{0, NEXT, PREV}, TransferInput(DIM_Y, 0, S_NEXT, S_PREV, shadow[0], D_NEXT, D_PREV, size[0], 1)},
        {{0, NEXT, NEXT}, TransferInput(DIM_Y, 0, S_NEXT, S_NEXT, shadow[0], D_NEXT, D_NEXT, size[0], 1)},
        {{PREV, PREV, PREV}, TransferInput(DIM_X, S_PREV, S_PREV, S_PREV, D_PREV, D_PREV, D_PREV, 1, 1)}, // corners XYZ
        {{PREV, PREV, NEXT}, TransferInput(DIM_X, S_PREV, S_PREV, S_NEXT, D_PREV, D_PREV, D_NEXT, 1, 1)},
        {{PREV, NEXT, PREV}, TransferInput(DIM_X, S_PREV, S_NEXT, S_PREV, D_PREV, D_NEXT, D_PREV, 1, 1)},
        {{PREV, NEXT, NEXT}, TransferInput(DIM_X, S_PREV, S_NEXT, S_NEXT, D_PREV, D_NEXT, D_NEXT, 1, 1)},
        {{NEXT, PREV, PREV}, TransferInput(DIM_X, S_NEXT, S_PREV, S_PREV, D_NEXT, D_PREV, D_PREV, 1, 1)},
        {{NEXT, PREV, NEXT}, TransferInput(DIM_X, S_NEXT, S_PREV, S_NEXT, D_NEXT, D_PREV, D_NEXT, 1, 1)},
        {{NEXT, NEXT, PREV}, TransferInput(DIM_X, S_NEXT, S_NEXT, S_PREV, D_NEXT, D_NEXT, D_PREV, 1, 1)},
        {{NEXT, NEXT, NEXT}, TransferInput(DIM_X, S_NEXT, S_NEXT, S_NEXT, D_NEXT, D_NEXT, D_NEXT, 1, 1)}
    };
    std::map<Index3, TransferInput> transfers;
    for (const auto &ni: neigh_indices.getNeighborIndices(params.neighbor_indices_type)) {
        auto it = all_transfers.find(ni);
        if (it != all_transfers.end()) {
            transfers.insert(*it);
        }
    }
    return transfers;
}

std::map<Index3, TransferInput> getBlockTransferInputs(const Index3 &block_index, const CellBlockParams &params,
                                                       const NeighborIndices &neigh_indices, const MeshParams &mesh_params) {
    auto transfer_inputs = getTransferInputs(params, neigh_indices);
    std::map<Index3, TransferInput> result;
    for (const auto &p: transfer_inputs) {
        Index3 neigh_index;
        for (int i = 0; i < 3; i++) {
            neigh_index[i] = (block_index[i] + mesh_params.num_blocks[i] + p.first[i]) % mesh_params.num_blocks[i]; // account periodic boundaries
        }
        result[neigh_index] = p.second;
    }
    return result;
}

struct Transfer {
    int dimension;
    int x, y, z;
    std::future<CellBlock::CellArray2> data;

    Transfer(int dim, int xx, int yy, int zz, std::future<CellBlock::CellArray2> &&dt) :
        dimension(dim), x(xx), y(yy), z(zz), data(std::move(dt)) {
    }

    void process(CellBlock &block) {
        ddl::Timer timer;
        auto cells = data.get();
        wait_time += timer.time();

        timer.reset();
        block.setSliceShadow(dimension, x, y, z, std::move(cells));
        slice_set_time += timer.time();
    }
};


void updateShadows(Mesh3D &mesh, const std::map<Index3, int> &block_distr,
                   const std::map<Index3, std::map<Index3, TransferInput>> &transfer_inputs) {
    static auto slice_get = [&](CellBlock &block, int dim, int x, int y, int z, int n1, int n2) -> CellBlock::CellArray2 {
        ddl::Timer timer;
        auto result = block.getSlice(dim, x, y, z, n1, n2);
        slice_get_time += timer.time();
        return result;
    };
    static auto get_slice = mesh.makeGetter<decltype(slice_get), CellBlock::CellArray2, int, int, int, int, int, int>(slice_get);

    std::vector<std::pair<CellBlock*, Transfer>> transfers;
    for (auto &p: mesh.localContents()) {
        auto index = p.first;
        auto *block = p.second.get();
        for (const auto &p: transfer_inputs.at(index)) {
            const auto &neigh_index = p.first;
            const auto &t_input = p.second;
            if (block_distr.find(neigh_index) != block_distr.end()) {
                auto future = get_slice(block_distr.at(neigh_index), neigh_index,
                                        t_input.dimension, t_input.src1, t_input.src2, t_input.src3,
                                        t_input.size1, t_input.size2);
                transfers.emplace_back(block, Transfer(t_input.dimension, t_input.dst1, t_input.dst2, t_input.dst3,
                                                       std::move(future)));
            }
        }
    }

    for (auto &p: transfers) {
        p.second.process(*p.first);
    }
}

int main(int argc, char **argv) {

    int size_x = intArg(1, 100, argc, argv);
    int size_y = intArg(2, 100, argc, argv);
    int size_z = intArg(3, 100, argc, argv);
    int num_of_particles = intArg(4, 1000, argc, argv);
    int num_of_iters = intArg(5, 100, argc, argv);
    //const bool scale_size = intArg(12, 0, argc, argv);
    //const bool scale_fragments = intArg(13, 0, argc, argv);
    //const bool debug = intArg(14, 0, argc, argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Comm cart_comm;

    int dims[3] = {0, 0, 0};
    MPI_Dims_create(size, 3, dims);
    int periods[3] = {0, 0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart_comm);

    /*// Size scaling mode: mesh size denotes mesh block size for each node.
    // Scale global mesh size appropriately.
    if (scale_size) {
        size_x *= nodes_by_x;
        size_y *= nodes_by_y;
        size_z *= nodes_by_z;
    }

    // Num of fragments scaling mode: num of fragments denotes num of fragments on each node.
    // Scale global num of fragments appropriately.
    if (scale_fragments) {
        num_of_fragments_x *= nodes_by_x;
        num_of_fragments_y *= nodes_by_y;
        num_of_fragments_z *= nodes_by_z;
    }*/

    MDParams md_params;
    md_params.epsilon = 5.0;
    md_params.sigma = 1.0;
    md_params.cutoff_radius = 2.5 * md_params.sigma;
    md_params.delta_t = 1e-5;

    md_params.load("conf.txt");

    const bool do_output = !md_params.output.empty();

    MeshParams mesh_params;
    mesh_params.mesh_size[0] = size_x;
    mesh_params.mesh_size[1] = size_y;
    mesh_params.mesh_size[2] = size_z;
    mesh_params.num_blocks[0] = dims[0];
    mesh_params.num_blocks[1] = dims[1];
    mesh_params.num_blocks[2] = dims[2];
    mesh_params.origin[0] = 0.0;
    mesh_params.origin[1] = 0.0;
    mesh_params.origin[2] = 0.0;
    mesh_params.step[0] = md_params.cutoff_radius;
    mesh_params.step[1] = md_params.cutoff_radius;
    mesh_params.step[2] = md_params.cutoff_radius;
    mesh_params.area_size[0] = mesh_params.step[0] * mesh_params.mesh_size[0];
    mesh_params.area_size[1] = mesh_params.step[1] * mesh_params.mesh_size[1];
    mesh_params.area_size[2] = mesh_params.step[2] * mesh_params.mesh_size[2];
    for (int i = 0; i < 3; i++) {
        mesh_params.periodic[i] = md_params.periodic[i];
    }

    MeshDecomosition mesh_decomp {{mesh_params.mesh_size[0], mesh_params.num_blocks[0]},
                                  {mesh_params.mesh_size[1], mesh_params.num_blocks[1]},
                                  {mesh_params.mesh_size[2], mesh_params.num_blocks[2]}};

    CellBlockParams block_params;
    for (size_t i = 0; i < 3; i++) {
        block_params.size[i] = mesh_decomp.blockSize(static_cast<int>(i), index[i]);
        block_params.shift[i] = mesh_decomp.blockShift(static_cast<int>(i), index[i]);
        block_params.shadow_start[i] = mesh_params.periodic[i] ? 1 : (index[i] > 0 ? 1 : 0);
        block_params.shadow_end[i] = mesh_params.periodic[i] ? 1 : (index[i] < mesh_params.num_blocks[i] - 1 ? 1 : 0);
        block_params.full_size[i] = block_params.size[i] + block_params.shadow_start[i] + block_params.shadow_end[i];
    }
    block_params.neighbor_indices_type = NeighborIndices::getType(
                block_params.shadow_start[0], block_params.shadow_end[0],
                block_params.shadow_start[1], block_params.shadow_end[1],
                block_params.shadow_start[2], block_params.shadow_end[2]);

    CellBlock block(block_params);
    //auto transfer_input = getBlockTransferInputs()
    NeighborIndices neigh_indices;

    std::map<Index3, std::map<Index3, TransferInput>> block_transfer_inputs;

    for (const auto &p: mesh.localContents()) {
        const auto index = p.first;
        const auto &params = p.second->getParams();
        block_transfer_inputs[index] = getBlockTransferInputs(index, params, neighbor_indices, mesh_params);
    }

    MPI_Barrier(cart_comm);

    initParticles([&](){
        return initParticles(num_of_particles, md_params, mesh_params);
    }, mesh, my_rank, mesh_params, mesh_decomp, distribution);

    if (do_output) {
        outputVTKNode(md_params.output, 0, my_rank, num_of_nodes, mesh);
    }

    env.getCommService()->getCommunicator()->barrier();

    double update_shadows_time = 0;
    double calc_forces_time = 0;
    double update_pos_time = 0;
    double update_vel_time = 0;
    double move_particles_time = 0;

    const auto t_start = MPI_Wtime();

    update_shadows_time += timed(updateShadows, mesh, distribution, block_transfer_inputs);

    calcForcesBlock(block, md_params, neigh_indices);

    MPI_Barrier(cart_comm);

    for (int iter = 1; iter <= num_of_iters; iter++) {
        updatePositionsBlock(block, md_params.delta_t);

        //env.getCommService()->getCommunicator()->barrier();

        move_particles_time += timed(moveParticles, mesh, mesh_params, mesh_decomp, distribution);

        //env.getCommService()->getCommunicator()->barrier();

        update_shadows_time += timed(updateShadows, mesh, distribution, block_transfer_inputs);

        //env.getCommService()->getCommunicator()->barrier();

        calcForcesBlock(block, md_params, neigh_indices);
        updateVelocitiesBlock(block, md_params.delta_t);

        if (do_output && (iter % md_params.output_step == 0)) {
            outputVTKNode(md_params.output, iter, rank, size, block);
        }
    }

    MPI_Barrier(cart_comm);

    const auto t_end = MPI_Wtime();

    const auto time = t_end - t_start;

    /*if (do_output) {
        outputVTKSeries(md_params.output, num_of_iters + 1, md_params);
    }*/

    std::ostringstream out;

    if (rank == 0) {
        out << "Mesh: " << mesh_params.mesh_size[0] << " x " << mesh_params.mesh_size[1] << " x " << mesh_params.mesh_size[2] <<
               ", fragments: " << mesh_params.num_blocks[0] << " x " << mesh_params.num_blocks[1] << " x " << mesh_params.num_blocks[2] <<
               ", grid: " << dims[0] << " x " << dims[1] << " x " << dims[2] <<
               ", particles: " << num_of_particles <<
               ", iters: " << num_of_iters <<
               ", time: " << time <<
               std::endl;
    }
    std::cout << out.str();

    return 0;
}
