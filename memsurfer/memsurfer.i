/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
Written by Harsh Bhatia (hbhatia@llnl.gov) and Peer-Timo Bremer (bremer5@llnl.gov)
LLNL-CODE-763493. All rights reserved.

This file is part of MemSurfer, Version 1.0.
Released under GNU General Public License 3.0.
For details, see https://github.com/LLNL/MemSurfer.
*/

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

%module memsurfer_cmod

%{
#define SWIG_FILE_WITH_INIT
#include "Types.hpp"
#include "PointSet.hpp"
#include "TriMesh.hpp"
#include "DensityKernels.hpp"
#include "DistanceKernels.hpp"
%}

%include "stdint.i"
%include "std_vector.i"
%include "std_string.i"

namespace std {
  %template(FloatVector) vector<float>;
  %template(IntVector) vector<int>;
}


%include "numpy.i"

%init %{
import_array();
%}

%apply (float* INPLACE_ARRAY1, int DIM1) {(float *_, int n)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *_, int n)};
%apply (float* INPLACE_ARRAY2, int DIM1, int DIM2) {(float *_, int n, int d)};
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *_, int n, int d)};
%apply (uint32_t* INPLACE_ARRAY2, int DIM1, int DIM2) {(uint32_t *_, int n, int d)};
%apply (int32_t* INPLACE_ARRAY2, int DIM1, int DIM2) {(int32_t *_, int n, int d)};

%apply (string key, float* INPLACE_ARRAY2, int DIM1, int DIM2) {(string key, float *_, int n, int d)};
%apply (string key, double* INPLACE_ARRAY2, int DIM1, int DIM2) {(string key, double *_, int n, int d)};

%include "Types.hpp"
%include "PointSet.hpp"
%include "TriMesh.hpp"
%include "DensityKernels.hpp"
%include "DistanceKernels.hpp"
