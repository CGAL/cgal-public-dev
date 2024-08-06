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

typedef CGAL::Surface_mesh_shortest_path_traits<K, Surface_mesh> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
typedef Traits::Barycentric_coordinates Barycentric_coordinates;

typedef boost::graph_traits<Surface_mesh> Graph_traits;
typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::face_iterator face_iterator;
typedef Graph_traits::vertex_descriptor vertex_descriptor;
typedef Graph_traits::face_descriptor face_descriptor;
typedef Graph_traits::halfedge_descriptor halfedge_descriptor;

typedef K::FT                                                           FT;
typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
//typedef PMP::Face_location<Surface_mesh, FT>                                    Face_location;
typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh>                AABB_face_graph_primitive;
typedef CGAL::AABB_traits_3<K, AABB_face_graph_primitive>               AABB_face_graph_traits;

typedef boost::graph_traits<Surface_mesh>::face_descriptor           face_descriptor;
typedef boost::graph_traits<Surface_mesh>::faces_size_type           faces_size_type;
typedef Surface_mesh::Property_map<face_descriptor, faces_size_type> FCCmap;
typedef CGAL::Face_filtered_graph<Surface_mesh>                      Filtered_graph;

typedef Surface_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;
typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Surface_mesh, CGAL::Heat_method_3::Direct> Heat_method_idt;

typedef PMP::Face_location<Surface_mesh, double>                      Face_location;
typedef PMP::Edge_location<Surface_mesh, double>                      Edge_location;

enum Distance_version { EUCLIDEAN_DISTANCE, GEODESIC_DISTANCE, HEAT_DISTANCE };

// A model of SurfacemeshShortestPathVisitor storing simplicies
// using boost::variant
struct Sequence_collector
{
  typedef boost::variant< vertex_descriptor,
                         std::pair<halfedge_descriptor,double>,
                         std::pair<face_descriptor, Barycentric_coordinates> > Simplex;
  std::vector< Simplex > sequence;
  void operator()(halfedge_descriptor he, double alpha)
  {
    sequence.push_back(std::make_pair(he, alpha));
  }
  void operator()(vertex_descriptor v)
  {
    sequence.push_back(v);
  }
  void operator()(face_descriptor f, Barycentric_coordinates alpha)
  {
    sequence.push_back(std::make_pair(f, alpha));
  }
};



double geodesicDistancePoints(Surface_mesh mesh, const Point source, const Point target, const CGAL::AABB_tree<AABB_face_graph_traits> &tree)
{

  Face_location query_location_source = PMP::locate_with_AABB_tree(source, tree, mesh);
  Face_location query_location_target = PMP::locate_with_AABB_tree(target, tree, mesh);

  // construct a shortest path query object and add a source point
   Surface_mesh_shortest_path shortest_paths(mesh);
  shortest_paths.add_source_point(query_location_source.first, {{query_location_source.second[0],query_location_source.second[1],query_location_source.second[2]}});

  // collect the sequence of simplicies crossed by the shortest path
  Sequence_collector sequence_collector;
  shortest_paths.shortest_path_sequence_to_source_points(query_location_target.first, {{query_location_target.second[0],query_location_target.second[1],query_location_target.second[2]}}, sequence_collector);

  std::vector<Point> points_in_path;
  shortest_paths.shortest_path_points_to_source_points(query_location_target.first, {{query_location_target.second[0],query_location_target.second[1],query_location_target.second[2]}}, std::back_inserter(points_in_path));
  /*Need to add the visitor from the geodesic distance example back in to print
  Print_visitor print_visitor(mesh);
  for (size_t i = 0; i < sequence_collector.sequence.size(); ++i)
    boost::apply_visitor(print_visitor, sequence_collector.sequence[i]);
  */
  double dist = 0;
  Point  pre_point = points_in_path[0];

  for (std::size_t i = 1; i < points_in_path.size(); ++i){
    dist = dist + sqrt(CGAL::squared_distance(points_in_path[i],pre_point));
    pre_point = points_in_path[i];
  }

  return dist;
}

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
double distancePoints(Surface_mesh mesh, const Point source, const Point target, const CGAL::AABB_tree<AABB_face_graph_traits> &tree, const Face_location start, const Face_location end)
{
  if constexpr (V==GEODESIC_DISTANCE)
    return geodesicDistancePoints(mesh, source, target, tree);
  if constexpr (V==EUCLIDEAN_DISTANCE)
    return euclideanDistancePoints(source, target);
  return 0;
}



std::vector<std::vector<std::vector<std::vector<Point>>>>
buildGrid(Surface_mesh /* mesh */, double minX, double maxX, double minY, double maxY, double minZ, double maxZ, double minDistance)
{
  // Define a four-dimensional vector
  std::vector<std::vector<std::vector<std::vector<Point>>>> four_d_grid;

  double cellSize = minDistance / sqrt(3.0);

  int gridX = (maxX - minX) / cellSize + 1;
  int gridY = (maxY - minY) / cellSize + 1;
  int gridZ = (maxZ - minZ) / cellSize + 1;

  for (int i = 0; i < gridX; ++i)
  {
    std::vector<std::vector<std::vector<Point>>> grid1;
    for (int j = 0; j < gridY; ++j)
    {
      std::vector<std::vector<Point>> grid2;
      for (int k = 0; k < gridZ; ++k)
      {
        std::vector<Point> grid3;
        grid2.push_back(grid3);
      }
      grid1.push_back(grid2);
    }
    four_d_grid.push_back(grid1);
  }

  return four_d_grid;
}

std::vector<std::vector<std::vector<std::vector<Point>>>>
addPointToGrid(Point point, std::vector<std::vector<std::vector<std::vector<Point>>>> grid, double minX, double minY, double minZ, double cellSize)
{
  int px =(1 / cellSize)*point.x() - minX*(1 / cellSize);
  int py =(1 / cellSize)*point.y() - minY*(1 / cellSize);
  int pz =(1 / cellSize)*point.z() - minZ*(1 / cellSize);

  grid[px][py][pz].push_back(point);

  return grid;
}


// Function to check if a point is within the minimum distance of existing points
template <Distance_version V>
bool
isFarEnoughFromExistingPoints(Surface_mesh mesh, Point point,
                              const std::vector<std::vector<std::vector<std::vector<Point>>>> grid,
                              double minX, double minY, double minZ, double minDistance, double cellSize/*, const CGAL::AABB_tree<AABB_face_graph_traits> &tree8 */)
{
  int gridX = (1 / cellSize)*point.x() - minX*(1 / cellSize);
  int gridY = (1 / cellSize)*point.y() - minY*(1 / cellSize);
  int gridZ = (1 / cellSize)*point.z() - minZ*(1 / cellSize);

  int minCheckX = fmax(0, gridX - 2);
  int maxCheckX = fmin(grid.size()-1 , gridX + 2);
  int minCheckY = fmax(0, gridY - 2);
  int maxCheckY = fmin(grid[0].size()-1 , gridY + 2);
  int minCheckZ = fmax(0, gridZ - 2);
  int maxCheckZ = fmin(grid[0][0].size()-1 , gridZ + 2);

  for (int i = minCheckX; i <= maxCheckX; ++i)
  {
    for (int j = minCheckY; j <= maxCheckY; ++j)
    {
      for(int k = minCheckZ; k <= maxCheckZ; ++k)
      {
        if(grid[i][j][k].size()==0){

        }
        else
        {
          for(Point p : grid[i][j][k])
          {
            if (euclideanDistancePoints(point, p) < minDistance)
                //distancePoints<V>(mesh, point, p, tree) < minDistance)
              return false;
          }
        }
      }
    }
  }
  return true;
}

// Generate  Sample
template <Distance_version V>
std::vector<Point> updatedPoissonDiskSampling(Surface_mesh mesh, std::size_t kMaxTries, double minDistance)
{
  double cellSize = minDistance / sqrt(3.0);

  std::array<Point, 8> obb_points;
   CGAL::oriented_bounding_box(mesh, obb_points,
                               CGAL::parameters::use_convex_hull(true));

  double maxX = obb_points[0].x();
  double minX = obb_points[0].x();
  double maxY = obb_points[0].y();
  double minY = obb_points[0].y();
  double maxZ = obb_points[0].z();
  double minZ = obb_points[0].z();
  for (const auto& point : obb_points)
  {
    if(point.x() > maxX)
        maxX = point.x();
    if(point.x() < minX)
        minX = point.x();
    if(point.y() > maxY)
        maxY = point.y();
    if(point.y() < minY)
        minY = point.y();
    if(point.z() > maxZ)
        maxZ = point.z();
    if(point.z() < minZ)
        minZ = point.z();
  }

  //CGAL::AABB_tree<AABB_face_graph_traits> tree;
  //PMP::build_AABB_tree(mesh, tree);

  std::vector<Point> points;
  std::queue<std::pair<Point, face_descriptor>> activePoints;

  CGAL::Random_points_in_triangle_mesh_3<Surface_mesh> g(mesh);
  Point c = *g;
    
  face_descriptor fc = g.last_item_picked();

  activePoints.push(std::make_pair(c, fc));
  points.push_back(c);

  Point currentPoint;
  face_descriptor currentFace;

  // Define a four-dimensional vector
  std::vector<std::vector<std::vector<std::vector<Point>>>> four_d_grid = buildGrid(mesh, minX, maxX, minY, maxY, minZ, maxZ, minDistance);

  four_d_grid = addPointToGrid(c, four_d_grid, minX, minY, minZ, cellSize);


  while (!activePoints.empty())
  {
      
    std::tie(currentPoint, currentFace) = activePoints.front();
    activePoints.pop();

    Face_location current_location = PMP::locate_in_face(currentPoint, currentFace, mesh);
   
    for (std::size_t i = 0; i < kMaxTries; ++i)
      {
        double angle = 2 * M_PI * (double)rand() / RAND_MAX;
        double distance = minDistance + minDistance * (double)rand() / RAND_MAX;
        K::Vector_2 dir(cos(angle),sin(angle));
        std::vector<Face_location> path = PMP::straightest_geodesic<K>(current_location, dir, distance, mesh);
        Point newPoint = PMP::construct_point(path.back(), mesh);

        if(isFarEnoughFromExistingPoints<V>(mesh, newPoint, four_d_grid, minX, minY, minZ, minDistance, cellSize))
        {
            activePoints.push(std::make_pair(newPoint,path.back().first));
            points.push_back(newPoint);
            four_d_grid = addPointToGrid(newPoint, four_d_grid, minX, minY, minZ, cellSize);
        }

    }

  }

  return points;
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
            Face_location old_location = PMP::locate_in_face(p, f, mesh);
            // if (distancePoints<V>(mesh, c, p, tree) < minDistance)
           // if (euclideanDistancePoints(c, p) < minDistance)
            if (geodesiceApproximation(c_location, old_location, mesh) < minDistance)
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
  std::vector<Point> points = updatedPoissonDiskSampling<EUCLIDEAN_DISTANCE>(mesh,kMaxTries,minDistance);
  timer.stop();
  std::cout << "Grid Euclidean done in " << timer.time() << "s.\n";
  std::ofstream out("sampling_euclidean.xyz");
  out << std::setprecision(17);
  std::copy(points.begin(), points.end(), std::ostream_iterator<Point>(out, "\n"));
  out.close();

  timer.reset();
  points.clear();
  timer.start();
  points = no_grid_Poisson_disk_sampling<EUCLIDEAN_DISTANCE>(mesh,kMaxTries,minDistance);
  timer.stop();
  std::cout << "No grid geodesic done in " << timer.time() << "s.\n";
  std::ofstream out1("sampling_no_grid.xyz");
  out1 << std::setprecision(17);
  std::copy(points.begin(), points.end(), std::ostream_iterator<Point>(out1, "\n"));
  out1.close();

    
  



/*

#if 1
  ///
  timer.reset();
  points.clear();
  timer.start();
  points = updatedPoissonDiskSampling<GEODESIC_DISTANCE>(mesh,kMaxTries,minDistance);
  timer.stop();
  std::cout << "Geodesic done in " << timer.time() << "s.\n";
  out.open("sampling_geodesic.xyz");
  out << std::setprecision(17);
  std::copy(points.begin(), points.end(), std::ostream_iterator<Point>(out, "\n"));
#endif
*/
  return 0;
}
