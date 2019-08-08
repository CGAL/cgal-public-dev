#ifndef ISR_BENCHMARK_DIST_UTILS_H
#define ISR_BENCHMARK_DIST_UTILS_H

// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

//Mesh
#include <CGAL/Surface_mesh.h>

//file includes
#include "isr_test_types.h"
#include "isr_test_util_reconstruction.h"

//boost
#include <boost/foreach.hpp>

//random
#include <CGAL/Random.h>

//AABB_tree
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

//AABB Tree
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;


// ----------------------------------------------------------------------------

class Measure_type
{
  public :

  Measure_type() {} ;

  std::string get_name_in_file() const {return(_name_in_file);}

  virtual double run(const Mesh &mesh, const PwnList &input_pwn) const {return 0.0;}

  protected :

  std::string _name_in_file;
};


class MeanDistPTM : public Measure_type
{
  public :

  MeanDistPTM() {} ;

  double run(const Mesh &mesh, const PwnList &input_pwn) const
  {
    //charging faces into AABB Tree
    Tree tree(faces(mesh).first, faces(mesh).second, mesh);

    //computation
    tree.accelerate_distance_queries();
    double sum;
    double sqd_dist;

    for (PwnList::const_iterator it = input_pwn.begin(); it != input_pwn.end(); ++it) {
      const Point& current_pt = it->first;
      sqd_dist = tree.squared_distance(current_pt);
      sum += CGAL::sqrt(sqd_dist);
    }
    double mean_dist = sum / (input_pwn.size());

    return( mean_dist );    
  }

  protected :

  std::string _name_in_file = "_mean_dist_ptm.dat" ;
};


#endif //ISR_BENCHMARK_DEST_UTILS_H