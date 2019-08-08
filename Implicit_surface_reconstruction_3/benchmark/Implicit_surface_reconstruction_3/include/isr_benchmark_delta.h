#ifndef ISR_BENCHMARK_DELTA_H
#define ISR_BENCHMARK_DELTA_H

// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

//Mesh
#include <CGAL/Surface_mesh.h>

//file includes
#include "isr_test_types.h"
#include "isr_test_util_bbox.h"

//boost
#include <boost/foreach.hpp>

//random
#include <CGAL/Random.h>


// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------

double compute_delta_position_noise(size_t lvl, double bbox)
{
  return (bbox * lvl / 1000);
}


#endif //ISR_BENCHMARK_DELTA_H