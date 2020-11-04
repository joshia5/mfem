//                                MFEM Example 1
//                              Omega_h Hello World
//
// Compile with: make ex1
//
// Sample runs:
//    ex1
//
// Note:         generates a 2x2 triangle mesh
//
// Description:  see note

#include "mfem.hpp"
#include <fstream>
#include <iostream>

#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_build.hpp>

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
