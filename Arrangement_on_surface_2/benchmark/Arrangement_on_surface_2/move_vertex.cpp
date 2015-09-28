#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/point_generators_2.h>

#include <iostream>
#include <fstream>
#include <boost/foreach.hpp>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Point_2 Point_2;
typedef CGAL::Random_points_in_square_2< Point_2 > Point_generator;
typedef CGAL::Arr_segment_traits_2<Kernel> Traits;
typedef CGAL::Arrangement_2<Traits> Arrangement_2;

int main()
{
  // testing random polygons
  CGAL::default_random = CGAL::Random(0);
  int nb_segments=1000;
  int nb_points=100;

  Point_generator ptgen(0.5);
  std::vector<Segment_2> segments(nb_segments);
  for (int i=0;i<nb_segments;++i)
    segments[i]=Segment_2(*ptgen++, *ptgen++);

  Arrangement_2 arr;
  insert(arr, segments.begin(), segments.end());

  // print the segments in a file
  std::ofstream output("input_segments.cgal");
  BOOST_FOREACH(const Segment_2& s, segments)
    output << "2 " << s[0] << " 0 " << s[1] << " 0\n";

  output.close();
  output.open("input_points.xyz");
  for(int i=0;i<nb_points; ++i)
  {
    Point_2 pt=*ptgen++;
    insert_point(arr, pt);
    output << pt << " 0\n";
  }
  output.close();

  std::cout << "Arrangement of " << nb_segments << " segments and " << nb_points << " points\n";
  std::cout << "arr.number_of_vertices() before " << arr.number_of_vertices() << "\n";

  std::vector<Arrangement_2::Vertex_handle> vertices;
  vertices.reserve(arr.number_of_vertices());
  for(Arrangement_2::Vertex_iterator it=arr.vertices_begin(),
                                     end=arr.vertices_end(); it!=end; ++it)
  {
    vertices.push_back( it );
  }

  CGAL::Timer time;

  time.start();
  int nb_topo_change=0;
  BOOST_FOREACH(Arrangement_2::Vertex_handle vh, vertices){
    bool change_in_topology=false;
    move_vertex(arr, vh, *ptgen++, change_in_topology);
    if (change_in_topology) ++nb_topo_change;
  }
  time.stop();
  std::cout << "arr.number_of_vertices() after random moves " << arr.number_of_vertices() << " - nb of topological change " << nb_topo_change << "\n";
  std::cout << "Time for random moves of vertices " << time.time() << "s" << std::endl;

  time.reset(); time.start();
  nb_topo_change=0;
  BOOST_FOREACH(Arrangement_2::Vertex_handle vh, vertices)
  {
    double xmove=CGAL::default_random.uniform_01<double>()/10000;
    double ymove=CGAL::default_random.uniform_01<double>()/10000;
    Point_2 new_pt(vh->point().x()+xmove, vh->point().y()+ymove);
    bool change_in_topology=false;
    move_vertex(arr, vh, new_pt, change_in_topology);
    if (change_in_topology) ++nb_topo_change;
  }
  time.stop();
  std::cout << "arr.number_of_vertices() after small moves " << arr.number_of_vertices() << " - nb of topological change " << nb_topo_change << "\n";
  std::cout << "Time for small moves of vertices " << time.time() << "s" << std::endl;
}
