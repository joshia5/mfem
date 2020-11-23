//                                MFEM Example 1
//                              Omega_h Hello World
//
// Compile with: make ex1
//
// Sample runs:
//    ex1
//
// coords_sizeote:         generates a 2x2 triangle mesh
//
// Description:  see note

#include "mfem.hpp"
//#include "unit_tests.hpp"
#include <fstream>
#include <iostream>

#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_build.hpp>

#include <Omega_h_for.hpp>
#include <Omega_h_element.hpp>

#ifndef MFEM_USE_OMEGAH
#error This example requires that MFEM is built with MFEM_USE_OMEGAH=YES
#endif

using namespace mfem;

int main(int argc, char *argv[])
{
  auto lib = Omega_h::Library();
  auto o_mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX,
    1., 1., 0, 2, 2, 0);
  Omega_h::vtk::write_parallel("box", &o_mesh);
  // things needed from omegaH mesh //
  int Dim = o_mesh.Omega_h::Mesh::dim();
  auto coords = o_mesh.Omega_h::Mesh::coords();
  const unsigned long int nverts = coords.size()/Dim;
  const unsigned long int nelems = o_mesh.Omega_h::Mesh::nelems();
  auto elem2vert = o_mesh.Omega_h::Mesh::ask_down(Dim, Omega_h::VERT);
  auto elem2vert_degree = Omega_h::element_degree(OMEGA_H_SIMPLEX,
    Dim, Omega_h::VERT);
  // to get the boundary and boundary elements, we will need to bring in the
  // ids of geom ents? or i think ther is an api which will give me the
  // classified elems.
  // can look at mark exposed sides and mark by exposure
  auto exposed_sides = mark_exposed_sides (o_mesh);
  auto ns = mesh->nents (mesh->dim() - 1); // num. of sides
  auto s2sc = mesh->ask_up (mesh->dim() - 1, mesh->dim()).a2ab;
  auto sc2c = mesh->ask_up (mesh->dim() - 1, mesh->dim()).ab2b;
  auto f = OMEGA_H_LAMBDA (LO s) {
    if ((s2sc[s + 1] - s2sc[s]) < 2) {
      //bdry elem = s2sc[s];
      ++NumOfBdrElements;
    }
  };
  parallel_for(ns, f, "count_bdrElems");
  Write<I8> boundary(NumOfBdrElements);// note the mfem boundary array is of
  // type Arrary<Element *>, so this boundary will need to be type cast
  // now get IDs of boundary elements
  int iter_bdrElems = 0;
  auto get_bdrElemId = OMEGA_H_LAMBDA (LO s) {
    if ((s2sc[s + 1] - s2sc[s]) < 2) {
      ++NumOfBdrElements;
      boundary[iter_bdrElems] = sc2c[s2sc[s]];// get the id of the side's
      // adjacent cell
    }
  };
  parallel_for(ns, get_bdrElemId, "get_bdrElemId");
  // after this get the verts of each doundary element using the ask_down
  // note that following the readPumiElement, 
  auto get_bdrVerts = OMEGA_H_LAMBDA (LO i) {
    boundary[iter_bdrElems] = sc2c[s2sc[s]];// get the id of the side's
  };
  parallel_for(NumOfBdrElements, get_bdrElemVerts, "get_bdrElemVerts");

  // now the ReadPumiElement, i.e. read elem2verts for BdrElements
  // note here an Element *el, which is ptr to mfem element will be created
  // and vertices will be assigned to it
  //note here the elements to vert connectivity will have to be created
  el = NewElement(classDim()); // ptr to mfem element
  int nv, *v;
  // Create element in MFEM
  nv = el->GetNVertices();
  v  = el->GetVertices();
  // Fill the connectivity
  for (int i = 0; i < nv; ++i) {
    v[i] = apf::getNumber(vert_num, Verts[i], 0, 0);
  }

  //v_num_loc is number of local vertices which apf creates. dont think its
  //needded from omegah mesh

  Device device("cuda");
  device.Print();
  const MemoryType d_mt = device.GetMemoryType();

  int coords_size = coords.size(); // Base vector size
  //int nSub = 5; // Sub vector size
  double *x = new double[coords_size]; // Allocate base vector data
  // Shallow copy into MFEM vector + read or write command to copy to the GPU
  Vector V(x, coords_size);
  V.ReadWrite();
  //V.Write();
  //V.Read();
  std::cout<<"V - Flags"<<std::endl;
  V.GetMemory().PrintFlags();

  // Shallow copy into MFEM vector + annotate as reference to a base
  // vector
  // Similar to doing double * Vs = &x[nSub]; // Pointer to sub-vector
  // data
  //Vector Vs;
  //Vs.MakeRef(V, nSub, (coords_size-nSub)); // MakeRef(vector, offset, size)

/*
  std::cout<<"Vs - Flags "<<std::endl;
  Vs.GetMemory().PrintFlags();
*/
  double *coords_mfem = V.ReadWrite(); // Pointer to device memory
  //double *d_x = V.Write(); // Pointer to device memory
  // double *d_x = V.Read(); // Pointer to device memory
  std::cout<<"Contents of V on the GPU"<<std::endl;

  auto fill = OMEGA_H_LAMBDA (Omega_h::LO i) {
    coords_mfem[i] = coords[i];
    printf("%.1f coords=%.1f\n", coords_mfem[i], coords[i]);
  };
  Omega_h::parallel_for(coords_size, fill);
  printf("\n");

  //printf("create mfem omgea mesh\n");
  //Omega_hMesh::Omega_hMesh(o_mesh, 0);
  //Mesh *mesh = Omega_hMesh(o_mesh, 0);

  return 0;
}
