#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/optimal_bounding_box.h>

#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <boost/property_map/property_map.hpp>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <CGAL/squared_distance_3.h>


#include <CGAL/Polygon_mesh_processing/Bsurf/locally_shortest_path.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Dynamic_property_map.h>

#include <CGAL/Real_timer.h>

#include <boost/variant.hpp>
#include <boost/lexical_cast.hpp>

#include <fstream>
#include <iostream>
#include <string>

#include <unordered_map>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <queue>



typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef K::Point_3                                                Point;
typedef K::Vector_3                                               Vector;
typedef CGAL::Surface_mesh<Point>                                 Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor      vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor        face_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

typedef boost::graph_traits<Surface_mesh> Graph_traits;
typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::face_iterator face_iterator;
typedef Graph_traits::vertex_descriptor vertex_descriptor;
typedef Graph_traits::face_descriptor face_descriptor;
typedef Graph_traits::halfedge_descriptor halfedge_descriptor;

typedef K::FT                                                           FT;
typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;


typedef boost::graph_traits<Surface_mesh>::face_descriptor           face_descriptor;
typedef boost::graph_traits<Surface_mesh>::faces_size_type           faces_size_type;
typedef Surface_mesh::Property_map<face_descriptor, faces_size_type> FCCmap;
typedef CGAL::Face_filtered_graph<Surface_mesh>                      Filtered_graph;

typedef Surface_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;
typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Surface_mesh, CGAL::Heat_method_3::Direct> Heat_method_idt;

typedef PMP::Face_location<Surface_mesh, double>                      Face_location;
typedef PMP::Edge_location<Surface_mesh, double>                      Edge_location;

enum Distance_version { EUCLIDEAN_DISTANCE, GEODESIC_DISTANCE };

double euclideanDistancePoints(const Point source, const Point target)
{
  return sqrt(CGAL::squared_distance(source,target));
}

double geodesiceApproximation(const Face_location source, const Face_location target, Surface_mesh mesh)
{
  PMP::Dual_geodesic_solver<double> solver;
  PMP::init_geodesic_dual_solver(solver, mesh);
  std::vector<Edge_location> edge_locations;
    
  CGAL::Polygon_mesh_processing::locally_shortest_path<double>(source, target, mesh, edge_locations, solver);

  return PMP::path_length<K>(edge_locations,source,target,mesh);
}

//function to switch between geodesic and Euclidean distance
template <Distance_version V>
double distancePoints(Surface_mesh mesh, const Point source, const Point target, const Face_location start, const Face_location end)
{
  if constexpr (V==GEODESIC_DISTANCE)
    return geodesiceApproximation(start, end, mesh);
  if constexpr (V==EUCLIDEAN_DISTANCE)
    return euclideanDistancePoints(source, target);
  return 0;
}


template <Distance_version V>
std::vector<face_descriptor> faces_in_sub_mesh(Point c, face_descriptor fc, Surface_mesh mesh, double minDistance)
{
    //Begin flooding
    face_descriptor fd = fc;

    std::vector<bool> selected(num_faces(mesh), false);
    std::vector<face_descriptor> selection;
    selected[fd] = true;
    selection.push_back(fd);

    auto do_queue_edge = [&](halfedge_descriptor h)
    {
      halfedge_descriptor hopp=opposite(h, mesh);
      if (is_border(hopp, mesh) || selected[face(hopp, mesh)]) return false;

      K::Segment_3 edge(mesh.point(source(h,mesh)), mesh.point(target(h,mesh)));

//      return (distancePoints<V>(mesh, mesh.point(source(h,mesh)), c, tree)< 3*minDistance ||
//              distancePoints<V>(mesh, mesh.point(target(h,mesh)), c, tree)< 3*minDistance);

      return (euclideanDistancePoints(mesh.point(source(h,mesh)), c)< 3*minDistance ||
              euclideanDistancePoints(mesh.point(target(h,mesh)), c)< 3*minDistance);
    };


    std::vector<halfedge_descriptor> queue;
    for (halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(fd, mesh), mesh))
      if (do_queue_edge(h))
        queue.push_back(opposite(h, mesh));

    while (!queue.empty())
    {
      halfedge_descriptor h = queue.back();
      face_descriptor f = face(h, mesh);
      queue.pop_back();
      if (!selected[f])
      {
        selected[f]=true;
        selection.push_back(f);
      }

      h=next(h, mesh);
      if (do_queue_edge(h)) queue.push_back(opposite(h, mesh));
      h=next(h, mesh);
      if (do_queue_edge(h)) queue.push_back(opposite(h, mesh));
    }

    /*
    std::cout << "center: " << c << "\n";
    for (face_descriptor f : selection)
      std::cout << f << " ";
    std::cout << "\n";
    */

    return selection;
}

template <Distance_version V>
bool
is_far_enough(Point c, Face_location c_location, Surface_mesh  mesh , double minDistance, std::vector<face_descriptor> selection, Surface_mesh::Property_map<face_descriptor,std::vector<Point>> face_points)
{
    for (face_descriptor f : selection)
    {
        for(Point p  : face_points[f])
        {
            Face_location p_location = PMP::locate_in_face(p, f, mesh);
            //Todo: Ask why  is distancePoints for Euclidean so much slower then just
            //calling euclideanDistancePoints?
            if (distancePoints<V>(mesh, c, p, c_location, p_location) < minDistance)
            //if (euclideanDistancePoints(c, p) < minDistance)
            //if (geodesiceApproximation(c_location, p_location, mesh) < minDistance)
                return false;
        }
    }
    return true;
}

template <Distance_version V>
std::vector<Point> no_grid_Poisson_disk_sampling(Surface_mesh mesh, std::size_t kMaxTries, double minDistance)
{
    std::vector<Point> points;
    std::queue<std::pair<Point, face_descriptor>> activePoints;

    CGAL::Random_points_in_triangle_mesh_3<Surface_mesh> g(mesh);
    Point c = *g;
    face_descriptor fc = g.last_item_picked();

    Surface_mesh::Property_map<face_descriptor,std::vector<Point>> face_points;
    face_points = mesh.add_property_map<face_descriptor, std::vector<Point>>("f:face_points").first;

    face_points[fc].push_back(c);
    activePoints.push(std::make_pair(c, fc));
    points.push_back(c);

    Point currentPoint;
    face_descriptor currentFace;
    while (!activePoints.empty())
    {
        std::tie(currentPoint, currentFace) = activePoints.front();
        activePoints.pop();

        std::vector<face_descriptor> selection = faces_in_sub_mesh<V>(currentPoint, currentFace, mesh, minDistance);
        Face_location current_location = PMP::locate_in_face(currentPoint, currentFace, mesh);
        int k = 0;
        while (k < kMaxTries)
        {
            double angle = 2 * M_PI * (double)rand() / RAND_MAX;
            double distance = minDistance + minDistance * (double)rand() / RAND_MAX;
            K::Vector_2 dir(cos(angle),sin(angle));
            std::vector<Face_location> path = PMP::straightest_geodesic<K>(current_location, dir, distance, mesh);
            
            Face_location newLocation = path.back();
            Point newPoint = PMP::construct_point(path.back(), mesh);

            if(is_far_enough<V>(newPoint, newLocation, mesh, minDistance, selection, face_points))
            {
                face_points[path.back().first].push_back(newPoint);
                activePoints.push(std::make_pair(newPoint,path.back().first));
                points.push_back(newPoint);
            }else
            {
                ++k;
            }

        }
    }

    return points;
}

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/eight.off");
  const double minDistance = (argc > 2) ? atof(argv[2]) : 0.05;
  Surface_mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  std::size_t kMaxTries = 30; // Number of attempts to find a point

  CGAL::Real_timer timer;
  timer.start();
  std::vector<Point> points = no_grid_Poisson_disk_sampling<EUCLIDEAN_DISTANCE>(mesh,kMaxTries,minDistance);
  timer.stop();
  std::cout << "Euclidean done in " << timer.time() << "s.\n";
  std::ofstream out("sampling_euclidean.xyz");
  out << std::setprecision(17);
  std::copy(points.begin(), points.end(), std::ostream_iterator<Point>(out, "\n"));
  out.close();



#if 1
  ///
  timer.reset();
  points.clear();
  timer.start();
  points = no_grid_Poisson_disk_sampling<GEODESIC_DISTANCE>(mesh,kMaxTries,minDistance);
  timer.stop();
  std::cout << "Geodesic done in " << timer.time() << "s.\n";
  out.open("sampling_geodesic.xyz");
  out << std::setprecision(17);
  std::copy(points.begin(), points.end(), std::ostream_iterator<Point>(out, "\n"));
#endif

  return 0;
}
