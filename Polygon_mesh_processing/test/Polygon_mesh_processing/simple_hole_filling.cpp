#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <vector>
#include <fstream>
#include <boost/tuple/tuple.hpp>

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Epic;


typedef typename  Epic::Point_3        Point_3;
typedef CGAL::Polyhedron_3<Epic>       Polyhedron;


void read_polyline_with_extra_points(
  const char* file_name, 
  std::vector<Point_3>& points,
  std::vector<Point_3>& extras)
{
  std::ifstream stream(file_name);
  if(!stream) { assert(false); }

  for(int i =0; i < 2; ++i) {
    int count;
    if(!(stream >> count)) { assert(false); }
    while(count-- > 0) {
      Point_3 p;
      if(!(stream >> p)) { assert(false); }
      i == 0 ? points.push_back(p) : extras.push_back(p);
    }
  }
}

void check_triangles(std::vector<Point_3>& points, std::vector<boost::tuple<int, int, int> >& tris) {
  if(points.size() - 3 != tris.size()) {
    std::cerr << "  Error: there should be n-2 triangles in generated patch." << std::endl;
    assert(false);
  }

  const int max_index = static_cast<int>(points.size())-1;
  for(std::vector<boost::tuple<int, int, int> >::iterator it = tris.begin(); it != tris.end(); ++it) {
    if(it->get<0>() == it->get<1>() ||
      it->get<0>() == it->get<2>() ||
      it->get<1>() == it->get<2>() ) 
    {
      std::cerr << "Error: indices of triangles should be all different." << std::endl;
      assert(false); 
    }  

    if(it->get<0>() >= max_index ||
      it->get<1>() >= max_index ||
      it->get<2>() >= max_index ) 
    {
      std::cerr << "  Error: max possible index check failed." << std::endl;
      assert(false);
    } 
  }
}


template<class HDS, class K>
class Polyhedron_builder : public CGAL::Modifier_base<HDS> {
  typedef typename K::Point_3 Point_3;
public:
  Polyhedron_builder(std::vector<boost::tuple<int, int, int> >* triangles, 
    std::vector<Point_3>* polyline) 
    : triangles(triangles), polyline(polyline) 
  { }

  void operator()(HDS& hds) {
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
    B.begin_surface(polyline->size() -1, triangles->size());

    for(typename std::vector<Point_3>::iterator it = polyline->begin();
      it != --polyline->end(); ++it) {
        B.add_vertex(*it);
    }

    for(typename std::vector<boost::tuple<int, int, int> >::iterator it = triangles->begin();
      it != triangles->end(); ++it) {
        B.begin_facet();
        B.add_vertex_to_facet(it->get<0>());
        B.add_vertex_to_facet(it->get<1>());
        B.add_vertex_to_facet(it->get<2>());
        B.end_facet();
    }

    B.end_surface();
  }

private:
  std::vector<boost::tuple<int, int, int> >* triangles;
  std::vector<Point_3>* polyline;
};

void check_constructed_polyhedron(const char* file_name,
  std::vector<boost::tuple<int, int, int> >* triangles, 
  std::vector<Point_3>* polyline,
  const bool save_poly) 
{
  Polyhedron poly;
  Polyhedron_builder<typename Polyhedron::HalfedgeDS,Epic> patch_builder(triangles, polyline);
  poly.delegate(patch_builder);

  if(!poly.is_valid()) {
    std::cerr << "  Error: constructed patch does not constitute a valid polyhedron." << std::endl;
    assert(false);
  }

  if (!save_poly)
    return;

  std::string out_file_name;
  out_file_name.append(file_name).append(".off");
  std::ofstream out(out_file_name.c_str());
  out << poly; out.close();
}


void kon_test(const char* file_name, bool use_DT, bool save_output) {
  std::vector<Point_3> points; // this will contain n and +1 repeated point
  std::vector<Point_3> extras;
  read_polyline_with_extra_points(file_name, points, extras);

  std::vector<boost::tuple<int, int, int> > tris;
  CGAL::Polygon_mesh_processing::triangulate_hole_polyline(
    points, extras, std::back_inserter(tris),
    CGAL::Polygon_mesh_processing::parameters::use_delaunay_triangulation(use_DT));

  check_triangles(points, tris);
  check_constructed_polyhedron(file_name, &tris, &points, save_output);

  std::cerr << "  Done!" << std::endl;
}



int main(int argc, char const *argv[])
{
	
  if(argc > 1)
  {
    const char* filename = argv[1];
    std::cout << "filename= " << filename << std::endl;
    bool use_DT = true;
    bool save_output = true;
    kon_test(filename, use_DT, save_output);
    
  }
  else
  {

  	std::cerr << "run what\n";
    
  }




	return 0;
}