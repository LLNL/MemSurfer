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

#include "TriMesh.hpp"

#ifdef VTK_AVAILABLE
#include "vtkPolyData.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkFloatArray.h"
#include "vtkCurvatures.h"
#endif

/// ----------------------------------------------------------------------------
//! compute gaussian and mean curvatures
/// ----------------------------------------------------------------------------

std::vector<TypeFunction> TriMesh::need_curvature(bool verbose) {

#ifndef VTK_AVAILABLE
    std::cerr << " ERROR: " << this->tag() << "::need_curvature - VTK not available! cannot compute curvatures!\n";
    return std::vector<TypeFunction>;
#else

    const size_t nverts = mVertices.size();

    if (mFields.find("curv_mean") == mFields.end()) {

        if (verbose) {
            std::cout << "   > " << this->tag() << "::need_curvature...";
            fflush(stdout);
        }

        std::vector<TypeFunction> curvature_mean(nverts);
        std::vector<TypeFunction> curvature_Gaussian(nverts);

        // create vtkpolydata object
        vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
        surface->Initialize();

        surface->SetPoints(vtkSmartPointer<vtkPoints>::New());
        surface->SetPolys(vtkSmartPointer<vtkCellArray>::New());
        surface->SetVerts(vtkSmartPointer<vtkCellArray>::New());

        for (auto it = mVertices.begin(); it != mVertices.end(); ++it) {
            const Vertex &v = *it;
            vtkIdType pid = surface->GetPoints()->InsertNextPoint(v[0], v[1], v[2]);
            surface->GetVerts()->InsertNextCell(1, &pid);
        }

        for (auto it = mFaces.begin(); it != mFaces.end(); ++it) {
            const Face &f = *it;
            vtkIdType cell[3] = {f[0], f[1], f[2]};
            surface->InsertNextCell(VTK_TRIANGLE,3,cell);
        }

        vtkSmartPointer<vtkCurvatures> mCurvaturesFilter =  vtkSmartPointer<vtkCurvatures>::New();
        mCurvaturesFilter->SetInputData(surface);
        mCurvaturesFilter->SetCurvatureTypeToMean();
        mCurvaturesFilter->Update();

        vtkSmartPointer<vtkCurvatures> gCurvaturesFilter =  vtkSmartPointer<vtkCurvatures>::New();
        gCurvaturesFilter->SetInputData(surface);
        gCurvaturesFilter->SetCurvatureTypeToGaussian();
        gCurvaturesFilter->Update();

        for(size_t i = 0; i < mVertices.size(); i++) {
            curvature_mean[i] = mCurvaturesFilter->GetOutput()->GetPointData()->GetScalars()->GetTuple(i)[0];
            curvature_Gaussian[i] = gCurvaturesFilter->GetOutput()->GetPointData()->GetScalars()->GetTuple(i)[0];
        }

        if(verbose)
            std::cout << " Done!\n";

        mFields["curv_mean"]  = curvature_mean;
        mFields["curv_gauss"] = curvature_Gaussian;
    }

    // return value!
    std::vector<TypeFunction> curvatures;
    curvatures.resize(2*nverts);

    const std::vector<TypeFunction> &curvature_mean = mFields["curv_mean"];
    const std::vector<TypeFunction> &curvature_Gaussian = mFields["curv_gauss"];

    for(size_t i = 0; i < nverts; i++) {
        curvatures[i]        = curvature_mean[i];
        curvatures[i+nverts] = curvature_Gaussian[i];
    }
    return curvatures;
#endif
}

/// ----------------------------------------------------------------------------
//! write as vtp file
/// ----------------------------------------------------------------------------
bool TriMesh::write_vtp(const std::string &fname) {

#ifndef VTK_AVAILABLE
    std::cerr << " ERROR: " << this->tag() << "::write_vtp - VTK not available!\n";
    return false;
#else
    std::cout << "   > " << this->tag() << "::write_vtp("<<fname<<")...";
    fflush(stdout);

    vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
    surface->Initialize();
    surface->SetPoints(vtkSmartPointer<vtkPoints>::New());
    surface->SetPolys(vtkSmartPointer<vtkCellArray>::New());

    // write vertices
    if (this->mDim == 2) {
        for (auto iter=mVertices.begin(); iter!=mVertices.end(); ++iter) {
            const Vertex &p = *iter;
            surface->GetPoints()->InsertNextPoint(p[0], p[1], 0.0);
        }
    }
    else {
        for (auto iter=mVertices.begin(); iter!=mVertices.end(); ++iter) {
            const Vertex &p = *iter;
            surface->GetPoints()->InsertNextPoint(p[0], p[1], p[2]);
        }
    }

    // write faces
    for (auto iter=mFaces.begin(); iter!=mFaces.end(); ++iter) {
        const Face &f = *iter;
        vtkIdType cell[3] = {f[0], f[1], f[2]};
        surface->InsertNextCell(VTK_TRIANGLE, 3, cell);
    }

    // write all the fields
    for (auto iter=mFields.begin(); iter!=mFields.end(); ++iter) {
        const std::string &name = iter->first;
        std::vector<TypeFunction> &data = iter->second;

        vtkSmartPointer<vtkFloatArray> field = vtkSmartPointer<vtkFloatArray>::New();
        field->SetNumberOfTuples(data.size());
        field->SetArray(data.data(), data.size(), 1);
        field->SetName(name.c_str());
        surface->GetPointData()->AddArray(field);
    }

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetInputData(surface);
    writer->SetFileName(fname.c_str());
    //writer->SetWriteToOutputString(true);
    writer->Write();

    std::cout << " Done! Wrote " << mVertices.size() << " vertices, "
                                 << mFaces.size() << " faces, and "
                                 << mFields.size() << " fields!\n";
    return true;
#endif
}

/// ----------------------------------------------------------------------------
/// ----------------------------------------------------------------------------
