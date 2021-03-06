#include "output.h"
#include "cell_block.h"
#include "md.h"

#include <fstream>

std::string pieceSource(const std::string &path, int iter, int this_node) {
    return path + "." + std::to_string(iter) + ".piece." + std::to_string(this_node) + ".vtp";
}

std::string parallelSource(const std::string &path, int iter) {
    return path + "." + std::to_string(iter) + ".pvtp";
}

std::string baseName(const std::string &filename) {
    auto base_name = filename;
    auto p = base_name.find_last_of("/");
    if (p == std::string::npos) {
        p = base_name.find_last_of("\\");
    }
    if (p != std::string::npos) {
        base_name = base_name.substr(p + 1);
    }
    return base_name;
}

void outputVTKPiece(const std::string &path, int iter, int this_node, const CellBlock &block) {
    const auto filename = pieceSource(path, iter, this_node);
    std::ofstream out(filename.c_str());
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    out << "<PolyData>\n";
    int num_of_particles = block.getNumOfParticles();
    out << "<Piece NumberOfPoints=\"" << num_of_particles << "\" " <<
           "NumberOfVerts=\"" << num_of_particles << "\" " <<
           "NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";

    out << "<Points>\n";
    out << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" name=\"Position\" format=\"ascii\">\n";
    block.forEachParticle([&out](const Particle &particle) {
        out << particle.pos.x << " " << particle.pos.y << " " << particle.pos.z << " ";
    });
    out << "\n</DataArray>\n";
    out << "</Points>\n";

    out << "<PointData Scalars=\"Mass\" Vectors=\"Velocity\">\n";
    out << "<DataArray type=\"Float64\" NumberOfComponents=\"1\" Name=\"Mass\" format=\"ascii\">\n";
    block.forEachParticle([&out](const Particle &particle) {
        out << particle.mass << " ";
    });
    out << "\n</DataArray>\n";
    out << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n";
    block.forEachParticle([&out](const Particle &particle) {
        out << particle.velocity.x << " " << particle.velocity.y << " " << particle.velocity.z << " ";
    });
    out << "\n</DataArray>\n";
    out << "</PointData>\n";

    out << "<Verts>\n";
    out << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 0; i < num_of_particles; i++) {
        out << i << " ";
    }
    out << "\n</DataArray>\n";
    out << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < num_of_particles; i++) {
        out << i + 1 << " ";
    }
    out << "\n</DataArray>\n";
    out << "</Verts>\n";

    out << "</Piece>\n";
    out << "</PolyData>\n";
    out << "</VTKFile>\n";
}

void outputVTKParallel(const std::string &path, int iter, int num_of_nodes) {
    const auto filename = parallelSource(path, iter);
    std::ofstream out(filename.c_str());
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"PPolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    out << "<PPolyData GhostLevel=\"0\">\n";
    out << "<PPoints>\n";
    out << "<PDataArray type=\"Float64\" NumberOfComponents=\"3\" name=\"Position\"/>\n";
    out << "</PPoints>\n";
    out << "<PPointData Scalars=\"Mass\" Vectors=\"Velocity\">\n";
    out << "<PDataArray type=\"Float64\" NumberOfComponents=\"1\" Name=\"Mass\"/>\n";
    out << "<PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\"/>\n";
    out << "</PPointData>\n";
    for (int i = 0; i < num_of_nodes; i++) {
        auto piece_source = baseName(pieceSource(path, iter, i));
        out << "<Piece Source=\"" << piece_source << "\"/>\n";
    }
    out << "</PPolyData>\n";
    out << "</VTKFile>\n";
}

void outputVTKNode(const std::string &path, int iter, int this_node, int num_of_nodes, const CellBlock &block) {
    outputVTKPiece(path, iter, this_node, block);
    if (this_node == 0) {
        outputVTKParallel(path, iter, num_of_nodes);
    }
}
