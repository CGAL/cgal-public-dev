#ifndef ISR_BENCHMARK_ARTEFACT_UTILS_H
#define ISR_BENCHMARK_ARTEFACT_UTILS_H

// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

//Mesh
#include <CGAL/Surface_mesh.h>

//file includes
#include "isr_test_types.h"
#include "isr_test_util_bbox.h"
#include "isr_benchmark_delta.h"

//boost
#include <boost/foreach.hpp>

//random
#include <CGAL/Random.h>


// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------

void position_noise_generator(const PwnList &input_pwnl, const size_t lvl, PwnList &modified_pwnl)
{
  double delta = compute_delta_position_noise(lvl, util_bb_diag(input_pwnl));

  BOOST_FOREACH(Point_with_normal pwn, input_pwnl) {
    Point p = pwn.first;
    Vector n = pwn.second;

    double new_x = p.x() + delta * CGAL::get_default_random().uniform_real<double>(-1.0,1.0);
    double new_y = p.y() + delta * CGAL::get_default_random().uniform_real<double>(-1.0,1.0);
    double new_z = p.z() + delta * CGAL::get_default_random().uniform_real<double>(-1.0,1.0);
    Point new_p(new_x, new_y, new_z);

    modified_pwnl.push_back(std::make_pair(new_p, n));
  }
}


#endif //ISR_BENCHMARK_ARTEFACT_UTILS_H