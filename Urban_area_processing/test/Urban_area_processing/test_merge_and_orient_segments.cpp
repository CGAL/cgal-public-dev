// STL includes.
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cassert>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Urban_area_processing/utils.h>

#include "include/Saver.h"
#include "include/Utilities.h"

namespace UAP = CGAL::Urban_area_processing;
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;
using Segment_map = CGAL::Identity_property_map<Segment_2>;
using Segments_2 = std::vector<Segment_2>;

using Saver = UAP::Saver<Kernel>;

// TODO:
// Improve this test by using a better set of initial segments.

int main(int argc, char *argv[]) {

  Saver saver;
  const std::string out_path = "/Users/monet/Documents/gf/urban-area-processing/logs/";
  std::cout << std::endl << "* Creating contours..." << std::endl;
  const FT scale = FT(1) / FT(5); // meters
  const FT min_length_2 = FT(1) / FT(2); // degrees

  Segment_map segment_map;
  Segments_2 segments = {
    Segment_2(Point_2(0, 0), Point_2(0.9, 0)),
    Segment_2(Point_2(1, 0), Point_2(1, 0.9)),
    Segment_2(Point_2(1, 1), Point_2(0.1, 1)),
    Segment_2(Point_2(0, 1), Point_2(0, 0.1))
  };
  saver.export_polylines(segments, 
    "/Users/monet/Documents/gf/urban-area-processing/logs/segments_initial");

  Kernel kernel;
  std::vector<Segments_2> contours;
  CGAL::Urban_area_processing::merge_and_orient_segments(
    kernel, segments, segment_map, 
    std::back_inserter(contours), 
    scale, min_length_2);

  std::cout << "Number of detected contours: " << contours.size() << std::endl;
  std::cout << std::endl;
  
  segments.clear();
  for (const auto& contour : contours)
    for (const auto& segment : contour)
      segments.push_back(segment);
  saver.export_polylines(segments, 
    "/Users/monet/Documents/gf/urban-area-processing/logs/segments_merged");

  return EXIT_SUCCESS;
}
