#pragma once

#include "mesh3d.h"

#include <string>

std::string pieceSource(const std::string &path, int iter, int this_node);
std::string parallelSource(const std::string &path, int iter);
void outputVTKPiece(const std::string &path, int iter, int this_node, const Mesh3D &mesh);
void outputVTKParallel(const std::string &path, int iter, int num_of_nodes);
void outputVTKSeries(const std::string &path, int num_of_iters, const MDParams &md_params);
void outputVTKNode(const std::string &path, int iter, int this_node, int num_of_nodes, const Mesh3D &mesh);
