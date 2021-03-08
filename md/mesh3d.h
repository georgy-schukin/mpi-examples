#pragma once

#include "didal.h"
#include "decomp/multi_block_decomposition.h"
#include "distr/mesh_block_distribution.h"
#include "index.h"
#include "md.h"

using Mesh3D = ddl::DistributedStorage<Index3, CellBlock>;
using MeshDecomosition = ddl::MultuBlockDecomposition<3>;
using MeshDistribution = ddl::MeshBlockDistribution<3>;
using BlockMap = std::map<Index3, int>;
