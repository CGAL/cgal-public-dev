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

//PMP
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

//dD Tree
#include <CGAL/point_generators_2.h>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

//AABB Tree
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

//dD_tree
typedef CGAL::Search_traits_3<Kernel> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree dD_Tree;

//dD_tree surface_mesh
typedef boost::graph_traits<Mesh>::vertex_descriptor dDT_Point;
typedef boost::graph_traits<Mesh>::vertices_size_type size_type;
typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type Vertex_point_pmap;
typedef CGAL::Search_traits_adapter<dDT_Point,Vertex_point_pmap,TreeTraits> SurfaceMeshTreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<SurfaceMeshTreeTraits> SurfaceMeshNeighbor_search;
typedef SurfaceMeshNeighbor_search::Tree SurfaceMeshdD_Tree;
typedef SurfaceMeshdD_Tree::Splitter Splitter;
typedef SurfaceMeshNeighbor_search::Distance Distance;

//Mesh index
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef Mesh::Halfedge_index halfedge_descriptor;
typedef Mesh::Vertex_index Vertex_index;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

// ----------------------------------------------------------------------------


class Measure_type
{
  public :

  Measure_type() {}

  std::string get_name_in_file() const {return(_name_in_file);}

  virtual double run(const Mesh &mesh, const PwnList &input_pwn) const {return 0.0;}

  protected :

  std::string _name_in_file;
};


class MeanDistPTM : public Measure_type
{
  public :

  MeanDistPTM() {_name_in_file = "_mean_dist_ptm.dat";}

  double run(const Mesh &mesh, const PwnList &input_pwn) const
  {
    //charging faces into AABB Tree
    Tree tree(faces(mesh).first, faces(mesh).second, mesh);

    //computation
    tree.accelerate_distance_queries();
    double sum = 0.0;

    for (PwnList::const_iterator it = input_pwn.begin(); it != input_pwn.end(); ++it) {
      const Point& current_pt = it->first;
      double sqd_dist = tree.squared_distance(current_pt);
      sum += CGAL::sqrt(sqd_dist);
    }
    double mean_dist = sum / (input_pwn.size());

    return( mean_dist );    
  }

};

class MeanDistMTP : public Measure_type
{
  public :

  MeanDistMTP() {_name_in_file = "_mean_dist_mtp.dat";}

  double run(const Mesh &mesh, const PwnList &input_pwn) const
  {
    //sampling mesh
    std::list<Point> sample_points;
    CGAL::Polygon_mesh_processing::sample_triangle_mesh(mesh,
                                                        std::back_inserter(sample_points),
                                                        4000);

    //putting input points into dD_Tree
    typedef typename PwnList::value_type PwnList_t;
    boost::function<Point(const PwnList_t&)> pwn_it_to_point_it = boost::bind(&PwnList_t::first, _1);
    dD_Tree tree(boost::make_transform_iterator(input_pwn.begin(), pwn_it_to_point_it), 
                  boost::make_transform_iterator(input_pwn.end(), pwn_it_to_point_it));

    //computation
    FT sum = 0;
    for (std::list<Point>::iterator it = sample_points.begin(); it != sample_points.end(); ++it) {
      Neighbor_search search(tree, *it, 1);
      sum += CGAL::sqrt(search.begin()->second);
    }
    double mean_dist = sum/(sample_points.size());

    return(mean_dist);
  }
};

class HausdorffPTM : public Measure_type
{
  public :

  HausdorffPTM() {_name_in_file = "_hausdorff_dist_ptm.dat";}

  double run(const Mesh &mesh, const PwnList &input_pwn) const
  {
    //charging faces into AABB Tree
    Tree tree(faces(mesh).first, faces(mesh).second, mesh);

    //computation
    tree.accelerate_distance_queries();
    FT sum;
    FT sqd_dist;
    FT max_sqd_dist = tree.squared_distance(input_pwn.begin()->first);

    for (PwnList::const_iterator it = input_pwn.begin(); it != input_pwn.end(); ++it) {
      const Point& current_pt = it->first;
      sqd_dist = tree.squared_distance(current_pt);
      if(sqd_dist > max_sqd_dist)
        max_sqd_dist = sqd_dist;
    }
    FT max_dist = CGAL::sqrt(max_sqd_dist);

    return(max_dist);
  }
};

class HausdorffMTP : public Measure_type
{
  public :

  HausdorffMTP() {_name_in_file = "_hausdorff_dist_mtp.dat";}

  double run(const Mesh &mesh, PwnList &input_pwn) const
  {
    typedef typename PwnList::value_type PwnList_t;
    boost::function<Point(const PwnList_t&)> pwn_it_to_point_it = boost::bind(&PwnList_t::first, _1);
    double max_dist = PMP::approximate_max_distance_to_point_set(
                                              mesh,
                                              CGAL::make_range( boost::make_transform_iterator(input_pwn.cbegin(), pwn_it_to_point_it),
                                                                 boost::make_transform_iterator(input_pwn.cend(), pwn_it_to_point_it)) ,
                                              4000 );

    return (max_dist);
  }
};

class MeanAngleDev : public Measure_type
{
  public :

  MeanAngleDev() {_name_in_file = "_mean_angle_dev.dat";}

  double run(Mesh &mesh, PwnList &input_pwn) const
  {
    double sum = 0.0;

    //compute each output mesh vertex's normal
    Mesh::Property_map<vertex_descriptor, Vector> vnormals_pm = 
      mesh.add_property_map<vertex_descriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;
    compute_area_weighted_vertex_normals(mesh, vnormals_pm);

    //putting mesh vertices into dD Tree
    Vertex_point_pmap vppmap = get(CGAL::vertex_point,mesh);
    SurfaceMeshdD_Tree tree(
              vertices(mesh).begin(),
              vertices(mesh).end(),
              Splitter(),
              SurfaceMeshTreeTraits(vppmap)
    );
    Distance tr_dist(vppmap);

    //for each input point, compute closest vertex of the mesh
    for(PwnList::const_iterator it = input_pwn.begin(); it != input_pwn.end(); ++it) {
      SurfaceMeshNeighbor_search search(tree, it->first, 1, 0, true, tr_dist);
      Vertex_index nearest_v = search.begin()->first;
      //compute deviation between input point normal and computed mesh vertex normal
      //add it up to the sum
      Vector in_normal = it->second;
      Vector out_normal = vnormals_pm[nearest_v];
      sum += std::acos( (CGAL::scalar_product(in_normal , out_normal)) / 
                          (CGAL::sqrt(in_normal.squared_length() * out_normal.squared_length())) );
    }
    return( sum/(input_pwn.size()) );
  }
};

class Measure_type_list
{
  public :

  Measure_type_list()
  {
    list.push_back(new MeanDistPTM());
    list.push_back(new MeanDistMTP());
    list.push_back(new HausdorffPTM());
    list.push_back(new HausdorffMTP());
    list.push_back(new MeanAngleDev());
  }

  ~Measure_type_list()
  {
    for (Measure_type* mt : list)
      delete mt;
  }


  std::list<Measure_type*> list;
};


#endif //ISR_BENCHMARK_DEST_UTILS_H