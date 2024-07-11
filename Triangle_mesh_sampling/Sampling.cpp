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
typedef PMP::Face_location<Surface_mesh, FT>                                    Face_location;

typedef boost::graph_traits<Surface_mesh>::face_descriptor           face_descriptor;
typedef boost::graph_traits<Surface_mesh>::faces_size_type           faces_size_type;
typedef Surface_mesh::Property_map<face_descriptor, faces_size_type> FCCmap;
typedef CGAL::Face_filtered_graph<Surface_mesh>                      Filtered_graph;

typedef Surface_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;
typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Surface_mesh, CGAL::Heat_method_3::Direct> Heat_method_idt;


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



double geodesicDistancePoints(Surface_mesh mesh, const Point source, const Point target)
{
  //Todo: compute tree first outside of function and reuse
  Face_location query_location_source = PMP::locate(source, mesh);
  Face_location query_location_target = PMP::locate(target, mesh);
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

  for (int i = 1; i < points_in_path.size(); ++i){
    dist = dist + sqrt(CGAL::squared_distance(points_in_path[i],pre_point));
    pre_point = points_in_path[i];
  }

  return dist;
}

double euclideanDistancePoints(Surface_mesh mesh, const Point source, const Point target)
{
  return sqrt(CGAL::squared_distance(source,target));
}

//function to switch between geodesic and Euclidean distance
template <Distance_version V>
double distancePoints(Surface_mesh mesh, const Point source, const Point target)
{
  if constexpr (V==GEODESIC_DISTANCE)
    return geodesicDistancePoints(mesh, source, target);
  if constexpr (V==EUCLIDEAN_DISTANCE)
    return euclideanDistancePoints(mesh, source, target);
  return 0;
}

// Generate sample of kMaxTries random points in the annulus around
// a given point.
template <Distance_version V>
std::vector<Point> randomPoints(Surface_mesh mesh, Point c, int kMaxTries, double minDistance)
{
  std::vector<Point> points;

  //Begin flooding
  //Todo: use tree only once
  Face_location c_location = PMP::locate(c, mesh);
  face_descriptor fd = c_location.first;

  std::vector<bool> selected(num_faces(mesh), false);
  std::vector<face_descriptor> selection;
  selected[fd] = true;
  selection.push_back(fd);

  auto do_queue_edge = [&](halfedge_descriptor h)
  {
    halfedge_descriptor hopp=opposite(h, mesh);
    if (is_border(hopp, mesh) || selected[face(hopp, mesh)]) return false;

    K::Segment_3 edge(mesh.point(source(h,mesh)), mesh.point(target(h,mesh)));
      //check if we want geodesic distance here
      //Todo: Update with distance to edge and use Euclidean distance
    return (distancePoints<V>(mesh, mesh.point(source(h,mesh)), c)< 2*minDistance ||
            distancePoints<V>(mesh, mesh.point(target(h,mesh)), c)< 2*minDistance);
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
  //create sub_mesh
  Surface_mesh::Property_map<face_descriptor,int> submap;
  submap = mesh.add_property_map<face_descriptor, int>("f:sub").first;
  for (face_descriptor f : selection)
      submap[f]=1;

  Surface_mesh sub_mesh;
  Filtered_graph ffg(mesh, 1, submap);
  copy_face_graph(ffg, sub_mesh);

  //std::cout << "Sub mesh:" << faces(sub_mesh).size() << std::endl;
  //Generate points in annulus
  //Todo: Consider isolated triangle problem
  //Todo: Polar geodesic sampling funtion here
  while(points.size() < kMaxTries)
  {
    CGAL::Random_points_in_triangle_mesh_3<Surface_mesh> b(sub_mesh);
    Point d = *b;
    double dis = distancePoints<V>(sub_mesh, d, c);
    if( dis > minDistance && dis < 2*minDistance)
      points.push_back(d);
  }

  return points;
}


std::vector<std::vector<std::vector<std::vector<Point>>>>
buildGrid(Surface_mesh mesh, double minX, double maxX, double minY, double maxY, double minZ, double maxZ, double minDistance)
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
                              double minX, double minY, double minZ, double minDistance, double cellSize)
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
            if (distancePoints<V>(mesh, point, p) < minDistance)
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
std::vector<Point> updatedPoissonDiskSampling(Surface_mesh mesh, int kMaxTries, double minDistance)
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

  std::vector<Point> points;
  std::queue<Point> activePoints;

  CGAL::Random_points_in_triangle_mesh_3<Surface_mesh> g(mesh);
  Point c = *g;
  activePoints.push(c);
  points.push_back(c);

  // Define a four-dimensional vector
  std::vector<std::vector<std::vector<std::vector<Point>>>> four_d_grid = buildGrid(mesh, minX, maxX, minY, maxY, minZ, maxZ, minDistance);

  four_d_grid = addPointToGrid(c, four_d_grid, minX, minY, minZ, cellSize);


  while (!activePoints.empty())
  {
    Point currentPoint = activePoints.front();
    activePoints.pop();

    std::vector<Point>  random_Points = randomPoints<V>(mesh, currentPoint, kMaxTries, minDistance);

    for (int i = 0; i < random_Points.size(); ++i){
      if(isFarEnoughFromExistingPoints<V>(mesh, random_Points[i], four_d_grid, minX, minY, minZ, minDistance, cellSize))
      {
        //std::cout << "I'm adding"<< std::endl;
        activePoints.push(random_Points[i]);
        points.push_back(random_Points[i]);
        four_d_grid = addPointToGrid(random_Points[i], four_d_grid, minX, minY, minZ, cellSize);
      }
    }
  }

  return points;
}


int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("../eight.off");
  const double minDistance = (argc > 2) ? atof(argv[2]) : 0.05;
  Surface_mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  int kMaxTries = 30; // Number of attempts to find a point

  CGAL::Real_timer timer;
  timer.start();
  std::vector<Point> points = updatedPoissonDiskSampling<EUCLIDEAN_DISTANCE>(mesh,kMaxTries,minDistance);
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
  points = updatedPoissonDiskSampling<GEODESIC_DISTANCE>(mesh,kMaxTries,minDistance);
  timer.stop();
  std::cout << "Geodesic done in " << timer.time() << "s.\n";
  out.open("sampling_geodesic.xyz");
  out << std::setprecision(17);
  std::copy(points.begin(), points.end(), std::ostream_iterator<Point>(out, "\n"));
#endif

  return 0;
}
