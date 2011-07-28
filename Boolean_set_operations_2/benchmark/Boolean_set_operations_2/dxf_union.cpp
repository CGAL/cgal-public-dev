/*! \file dxf_union.cpp
 * Computing the union of a set of circular polygons read from a DXF file.
 */

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/IO/Dxf_bsop_reader.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <time.h>

//to enable user to enable or disable OpenMP only for Boolean_set_operations_2 package
#ifndef CGAL_BOOLEAN_SET_OPERATIONS_2_DONT_USE_OPENMP
#define CGAL_BOOLEAN_SET_OPERATIONS_2_USE_OPENMP 
#endif

//add OpenMP header
#if defined _OPENMP && defined CGAL_BOOLEAN_SET_OPERATIONS_2_USE_OPENMP
#include<omp.h>
#endif

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Gps_circle_segment_traits_2<Kernel>    Traits_2;
typedef Traits_2::Polygon_2                          Circ_polygon_2;
typedef Traits_2::Polygon_with_holes_2               Circ_polygon_with_holes_2;
typedef std::vector<Circ_polygon_2>                  Polygons_vec;
typedef std::vector<Circ_polygon_with_holes_2>       Polygons_with_holes_vec;
typedef CGAL::General_polygon_set_2<Traits_2>        General_polygon_set_2;

static const int DEFAULT_GROUP_SIZE = 5;

// The command line should be:
//   ex_dxf_union [DXF file] [simplify] [group size]
int main (int argc, char* argv[])
{
  // Open the input DXF file.
  const char* filename = (argc >= 2) ? argv[1] : "test.dxf";
  std::ifstream input_file (filename);
  if (! input_file.is_open())
  {
    std::cerr << "Failed to open the " << filename <<std::endl;
    return -1;
  }

  // Read the extra flags.
  bool    simplify = true;
  int     group_size = DEFAULT_GROUP_SIZE;

  if (argc >= 3)
  {
    simplify = (std::atoi (argv[2]) != 0);

    if (argc >= 4)
      group_size = std::atoi (argv[3]);

    if (group_size <= 0)
    {
      std::cerr << "Illegal group size: " << group_size << "." << std::endl;
      return (1);
    }
  }

  // Read the circular polygons from the DXF file.
  General_polygon_set_2          gps;
  Polygons_vec                   pgns;
  Polygons_with_holes_vec        pgns_with_holes;
  CGAL::Dxf_bsop_reader<Kernel>  reader;
  CGAL::Timer                    t_read;

  std::cout << "Reading <" << filename << "> ... " << std::flush;
  t_read.start();

  reader (input_file,
          std::back_inserter(pgns),
          std::back_inserter(pgns_with_holes),
          simplify);

  t_read.stop();
  std::cout << "Done! (" << t_read.time() << " seconds)." << std::endl;
  std::cout << std::distance (pgns.begin(), pgns.end()) << " polygons, "
            << std::distance (pgns_with_holes.begin(), pgns_with_holes.end())
            << " polygons with holes." << std::endl;

  input_file.close();

  // Compute their union.
  CGAL::Timer                    t_union;
  time_t                         s,e;

  std::cout << "Computing the union ... " << std::flush;

  t_union.start();
  s = time(NULL);

  #if defined _OPENMP && defined CGAL_BOOLEAN_SET_OPERATIONS_2_USE_OPENMP
  double start,end;
  start = omp_get_wtime(); 
  #endif     
	
  gps.join (pgns.begin(), pgns.end(),
            pgns_with_holes.begin(), pgns_with_holes.end(),
            group_size);
 
  #if defined _OPENMP && defined CGAL_BOOLEAN_SET_OPERATIONS_2_USE_OPENMP
  end = omp_get_wtime();
  #endif

  t_union.stop();
  e = time(NULL);

  #if defined _OPENMP && defined CGAL_BOOLEAN_SET_OPERATIONS_2_USE_OPENMP
  std::cout << "Done! time using omp_get_wtime() (" << end-start << " seconds)." << std::endl;
  #endif 
  
  std::cout << "Done! time using CGAL::Timer (using measure of clock cycles) (" << t_union.time() << " seconds).\n";

  std::cout << "Done! time using time() (" << difftime(e,s) << " seconds)." << std::endl;

  std::cout << "The result:"
            << "  |V| = " << gps.arrangement().number_of_vertices()
            << "  |E| = " << gps.arrangement().number_of_edges()
            << "  |F| = " << gps.arrangement().number_of_faces() << std::endl;

  return (0);
}
