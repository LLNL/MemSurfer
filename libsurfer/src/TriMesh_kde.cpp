/**
Copyright (c) 2019, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
Written by Harsh Bhatia (hbhatia@llnl.gov) and Peer-Timo Bremer (bremer5@llnl.gov)
LLNL-CODE-763493. All rights reserved.

This file is part of MemSurfer, Version 1.0.
Released under GNU General Public License 3.0.
For details, see https://github.com/LLNL/MemSurfer.
*/

/// ----------------------------------------------------------------------------
/// ----------------------------------------------------------------------------

#include <map>
#include <sstream>
#include <stdexcept>

#include "TriMesh.hpp"
#include "DensityEstimation.hpp"

/// -----------------------------------------------------------------------------
//! density estimation (core function)
/// -----------------------------------------------------------------------------

void kde_2d(const std::vector<Vertex> &vertices, const std::vector<TypeIndexI> &ids,
              const DensityKernel& k, const DistanceKernel &dist,
              std::vector<TypeFunction> &density) {

    size_t nids = ids.size();
    size_t nverts = vertices.size();
    density.resize(nverts, 0);

    if (ids.empty()) {  // compute for all ids
        TypeFunction tmp;
        for (TypeIndex j=0;   j<nverts; j++) {
        for (TypeIndex i=j+1; i<nverts; i++) {

            tmp = k(dist(vertices[i][0], vertices[i][1], vertices[j][0], vertices[j][1]));

            density[j] += tmp;
            density[i] += tmp;
        }}
    }
    else {              // compute for selected ids
        for (TypeIndex j=0; j<nverts; j++) {
        for (TypeIndex i=0; i<nids;   i++) {

            density[j] += k(dist(vertices[ids[i]][0], vertices[ids[i]][1], vertices[j][0], vertices[j][1]));
        }}
    }
}

void kde_3d(const std::vector<Vertex> &vertices, const std::vector<TypeIndexI> &ids,
              const DensityKernel& k, const DistanceKernel &dist,
              std::vector<TypeFunction> &density) {

    size_t nids = ids.size();
    size_t nverts = vertices.size();
    density.resize(nverts, 0);

    if (ids.empty()) {  // compute for all ids
        TypeFunction tmp;
        for (TypeIndex j=0;   j<nverts; j++) {
        for (TypeIndex i=j+1; i<nverts; i++) {

            tmp = k(dist(vertices[i][0], vertices[i][1], vertices[i][2],
                         vertices[j][0], vertices[j][1], vertices[j][2]));

            density[j] += tmp;
            density[i] += tmp;
        }}
    }
    else {              // compute for selected ids
        for (TypeIndex j=0; j<nverts; j++) {
        for (TypeIndex i=0; i<nids;   i++) {

            density[j] += k(dist(vertices[ids[i]][0], vertices[ids[i]][1], vertices[ids[i]][2],
                                 vertices[j][0], vertices[j][1], vertices[j][2]));
        }}
    }
}

/// -----------------------------------------------------------------------------
//! density estimation for nonperiodic mesh using all vertices
/// -----------------------------------------------------------------------------

const std::vector<TypeFunction>& TriMesh::kde(const DensityKernel& k, const std::string &name,
                                              bool verbose) {

    // this name already exists in the map!
    if (mFields.find(name) != mFields.end()) {
        if (verbose)
            std::cout << " TriMesh::kde("<<name<<") got an already used name. Returning existing field!\n";
        return mFields.at(name);
    }

    if (verbose){
        std::cout << "   > TriMesh::kde("<<name<<")...";
        fflush(stdout);
    }

    DistanceSquared dist;
    kde_2d(mVertices, std::vector<TypeIndexI>(), k, dist, mFields[name]);

    if (verbose){
        printf(" Done!\n");
    }
    return mFields.at(name);
}

/// -----------------------------------------------------------------------------
//! density estimation for nonperiodic mesh using a subset of vertices
/// -----------------------------------------------------------------------------

const std::vector<TypeFunction>& TriMesh::kde(const DensityKernel& k, const std::string &name,
                                              const std::vector<TypeIndexI> &ids, bool verbose) {

    // this name already exists in the map!
    if (mFields.find(name) != mFields.end()) {
        if (verbose){
            std::cout << " " << this->mName << "::kde("<<name<<") got an already used name. Returning existing field!\n";
        }
        return mFields.at(name);
    }

    if (verbose){
        std::cout << "   > " << this->mName << "::kde(<" << ids.size() << ">, "<<name<<")...";
        fflush(stdout);
    }

    DistanceSquared dist;
    kde_2d(mVertices, ids, k, dist, mFields[name]);

    if(verbose){
        printf(" Done!\n");
    }
    return mFields.at(name);
}

/// -----------------------------------------------------------------------------
//! density estimation for periodic mesh using all of vertices
/// -----------------------------------------------------------------------------

const std::vector<TypeFunction>& TriMeshPeriodic::kde(const DensityKernel& k, const std::string &name,
                                                      bool verbose) {

    // this name already exists in the map!
    if (mFields.find(name) != mFields.end()) {
        if (verbose)
            std::cout << " TriMesh::kde("<<name<<") got an already used name. Returning existing field!\n";
        return mFields.at(name);
    }

    if (verbose){
        std::cout << "   > TriMesh::kde("<<name<<")...";
        fflush(stdout);
    }

    DistancePeriodicXYSquared dist(mBox0, mBox1);
    kde_2d(mVertices, std::vector<TypeIndexI>(), k, dist, mFields[name]);

    if (verbose){
        printf(" Done!\n");
    }
    return mFields.at(name);
}

/// -----------------------------------------------------------------------------
//! density estimation for periodic mesh using a subset of vertices
/// -----------------------------------------------------------------------------

const std::vector<TypeFunction>& TriMeshPeriodic::kde(const DensityKernel& k, const std::string &name,
                                                      const std::vector<TypeIndexI> &ids, bool verbose) {

    // this name already exists in the map!
    if (mFields.find(name) != mFields.end()) {
        if (verbose){
            std::cout << " " << this->mName << "::kde("<<name<<") got an already used name. Returning existing field!\n";
        }
        return mFields.at(name);
    }

    if (verbose){
        std::cout << "   > " << this->mName << "::kde(<" << ids.size() << ">, "<<name<<")...";
        fflush(stdout);
    }

    DistancePeriodicXYSquared dist(mBox0, mBox1);
    kde_2d(mVertices, ids, k, dist, mFields[name]);

    if(verbose){
        printf(" Done!\n");
    }
    return mFields.at(name);
}

/// -----------------------------------------------------------------------------
/// -----------------------------------------------------------------------------
