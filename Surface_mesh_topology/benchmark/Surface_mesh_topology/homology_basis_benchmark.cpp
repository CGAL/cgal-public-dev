// Copyright (c) 2026 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Raphael Costil
//
// Benchmarks Minimal_length_homology_basis::compute_basis() across a set
// of meshes, repeating each one a fixed number of times, and writes one
// CSV row per run to the output file: label,expected_genus,computed_genus,V,F,run_index,time_seconds.
//
// Reads its mesh list from a manifest file (see
// prepare_homology_basis_meshes.py), not from argv directly -- one line
// per mesh, "label,genus,path", with a header row.
//
// Usage: homology_basis_benchmark <manifest.csv> <repeats> <output.csv>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Timer.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using K = CGAL::Simple_cartesian<double>;
using Mesh = CGAL::Surface_mesh<K::Point_3>;

struct Mesh_entry
{
  std::string label;
  int expected_genus;
  std::string path;
};

std::vector<Mesh_entry> read_manifest(const std::string& manifest_path)
{
  std::vector<Mesh_entry> entries;
  std::ifstream in(manifest_path);
  if (!in)
  {
    std::cerr << "ERROR: cannot open manifest " << manifest_path << std::endl;
    return entries;
  }

  std::string line;
  std::getline(in, line); // header
  while (std::getline(in, line))
  {
    if (!line.empty())
    {
      std::istringstream iss(line);
      std::string label, genus_str, path;
      std::getline(iss, label, ',');
      std::getline(iss, genus_str, ',');
      std::getline(iss, path, ',');
      entries.push_back({label, std::stoi(genus_str), path});
    }
  }
  return entries;
}

int main(int argc, char** argv)
{
  if (argc != 4)
  {
    std::cerr << "Usage: " << argv[0] << " <manifest.csv> <repeats> <output.csv>\n";
    return EXIT_FAILURE;
  }

  std::vector<Mesh_entry> entries = read_manifest(argv[1]);
  int repeats = std::stoi(argv[2]);

  std::ofstream out(argv[3]);
  if (!out)
  {
    std::cerr << "ERROR: cannot open output file " << argv[3] << std::endl;
    return EXIT_FAILURE;
  }
  out << "label,expected_genus,computed_genus,V,F,run_index,time_seconds\n";

  for (const Mesh_entry& entry : entries)
  {
    Mesh mesh;
    if (!CGAL::IO::read_polygon_mesh(entry.path, mesh))
    {
      std::cerr << "ERROR reading file " << entry.path << std::endl;
    }
    else
    {
      for (int run=0; run<repeats; ++run)
      {
        CGAL::Surface_mesh_topology::Curves_on_surface_topology<Mesh> cst(mesh);
        CGAL::Surface_mesh_topology::Euclidean_length_weight_functor<Mesh> wf(mesh);

        CGAL::Timer t;
        t.start();
        auto basis = cst.compute_minimal_homology_basis(wf);
        t.stop();

        int computed_genus = static_cast<int>(basis.size()) / 2;

        out << entry.label << ','
            << entry.expected_genus << ','
            << computed_genus << ','
            << mesh.number_of_vertices() << ','
            << mesh.number_of_faces() << ','
            << run << ','
            << t.time() << '\n';
        out.flush();

        std::cout << "  [" << entry.label << "] run " << run << "/" << repeats
                   << ": " << t.time() << " s" << std::endl;
      }
    }
  }

  return EXIT_SUCCESS;
}
