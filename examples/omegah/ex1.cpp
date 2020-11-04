//                                MFEM Example 1
//                              PUMI Modification
//
// Compile with: make ex1
//
// Sample runs:
//    ex1 -m ../../data/pumi/serial/Kova.smb -p ../../data/pumi/geom/Kova.dmg
//
// Note:         Example models + meshes for the PUMI examples can be downloaded
//               from github.com/mfem/data/pumi. After downloading we recommend
//               creating a symbolic link to the above directory in ../../data.
//
// Description:  This example code demonstrates the use of MFEM to define a
//               simple finite element discretization of the Laplace problem
//               -Delta u = 1 with homogeneous Dirichlet boundary conditions.
//               Specifically, we discretize using a FE space of the specified
//               order, or if order < 1 using an isoparametric/isogeometric
//               space (i.e. quadratic for quadratic curvilinear mesh, NURBS for
//               NURBS mesh, etc.)
//
//               The example highlights the use of mesh refinement, finite
//               element grid functions, as well as linear and bilinear forms
//               corresponding to the left-hand side and right-hand side of the
//               discrete linear system. We also cover the explicit elimination
//               of essential boundary conditions, static condensation, and the
//               optional connection to the GLVis tool for visualization.
//
//               This PUMI modification demonstrates how PUMI's API can be used
//               to load a PUMI mesh classified on a geometric model and then
//               convert it to the MFEM mesh format. The inputs are a Parasolid
//               model, "*.xmt_txt" and a SCOREC mesh "*.smb". The option "-o"
//               is used for the Finite Element order and "-go" is used for the
//               geometry order. Note that they can be used independently, i.e.
//               "-o 8 -go 3" solves for 8th order FE on a third order geometry.
//
// NOTE:         Model/Mesh files for this example are in the (large) data file
//               repository of MFEM here https://github.com/mfem/data under the
//               folder named "pumi", which consists of the following sub-folders:
//               a) geom -->  model files
//               b) parallel --> parallel pumi mesh files
//               c) serial --> serial pumi mesh files

#include "mfem.hpp"
#include <fstream>
#include <iostream>

#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_comm.hpp>

#ifndef MFEM_USE_OMEGAH
#error This example requires that MFEM is built with MFEM_USE_OMEGAH=YES
#endif

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
   auto lib = Omega_h::Library();
   auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX,
       1., 1., 0, 2, 2, 0);
   Omega_h::vtk::write_parallel("box", &mesh);
   return 0;
}
