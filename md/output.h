#pragma once

class CellBlock;
class MDParams;

#include <string>

std::string pieceSource(const std::string &path, int iter, int this_node);
std::string parallelSource(const std::string &path, int iter);
void outputVTKPiece(const std::string &path, int iter, int this_node, const CellBlock &block);
void outputVTKParallel(const std::string &path, int iter, int num_of_nodes);
void outputVTKNode(const std::string &path, int iter, int this_node, int num_of_nodes, const CellBlock &block);
