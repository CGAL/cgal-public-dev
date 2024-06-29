#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/Polygon_mesh_processing/locate.h>


#include <boost/variant.hpp>
#include <boost/lexical_cast.hpp>

#include <fstream>
#include <iostream>
#include <string>

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
// A visitor to print what a variant contains using boost::apply_visitor
struct Print_visitor
  : public boost::static_visitor<>
{
  int i;
  Surface_mesh& g;
  Print_visitor(Surface_mesh& g) :i(-1), g(g) {}
  void operator()(vertex_descriptor v)
  {
    std::cout << "#" << ++i << " Vertex: " << get(boost::vertex_index, g)[v];
    std::cout << " Position: " << Surface_mesh_shortest_path::point(v, g) << "\n";
  }
  void operator()(const std::pair<halfedge_descriptor,double>& h_a)
  {
    std::cout << "#" << ++i << " Edge: " << get(CGAL::halfedge_index, g)[h_a.first] << " , ("
                                           << 1.0 - h_a.second << " , "
                                           << h_a.second << ")";
    std::cout << " Position: " << Surface_mesh_shortest_path::point(h_a.first, h_a.second, g) << "\n";
  }
  void operator()(const std::pair<face_descriptor, Barycentric_coordinates>& f_bc)
  {
    std::cout << "#" << ++i << " Face: " << get(CGAL::face_index, g)[f_bc.first] << " , ("
                                           << f_bc.second[0] << " , "
                                           << f_bc.second[1] << " , "
                                           << f_bc.second[2] << ")";
    std::cout << " Position: " << Surface_mesh_shortest_path::point(f_bc.first, f_bc.second, g) << "\n";
  }
};



int kMaxTries = 30; // Number of attempts to find a point
double kMinDistance = 20.0; // Minimum distance between points


// Function to check if a point is in a rectangle
bool isInRectangle(double x, double y, double rectWidth, double rectHeight) {
    return x >= 0 && x <= rectWidth && y >= 0 && y <= rectHeight;
}

// Function to check if a point is in a triangle
bool isInTriangle(double x, double y, double x1, double x2, double y2) {
    return y >= 0 && y <= (y2/x2)*x && y <= (y2/(x2-x1))*(x-x1);
}




/*
// Generate  Sampling Rectangle
std::vector<Point> generatePoissonDiskSampling(double x1, double x2, double y2, double minDistance) {

    double width = x1;
    double height = abs(y2);
    srand(time(NULL));
    double cellSize = minDistance / sqrt(2.0);

    int gridWidth = (int)ceil(width / cellSize);
    int gridHeight = (int)ceil(height / cellSize);

    std::vector<std::vector<bool>> grid(gridWidth, std::vector<bool>(gridHeight, false));
    std::vector<Point> points;
    std::queue<Point> activePoints;

    // Generate first point
    //double startX = width * (double)rand() / RAND_MAX;
    //double startY = height * (double)rand() / RAND_MAX;
    Point random = randomInTriangle(x1,x2,y2);
    std::cout << "Random Point: (" << random.x() << ", " << random.y() <<  ")" << std::endl;
    double startX = random.x();
    double startY = random.y();
    activePoints.push(Point(startX, startY, 0));
    points.push_back(Point(startX, startY, 0));
    grid[(int)(startX / cellSize)][(int)(startY / cellSize)] = true;

    while (!activePoints.empty()) {
        Point currentPoint = activePoints.front();
        activePoints.pop();

        for (int i = 0; i < kMaxTries; ++i) {
            double angle = 2 * M_PI * (double)rand() / RAND_MAX;
            double distance = minDistance + minDistance * (double)rand() / RAND_MAX;
            double newX = currentPoint.x() + distance * cos(angle);
            double newY = currentPoint.y() + distance * sin(angle);
           // std::cout << "Test points: (" << newX << ", " << newY <<  ")" << std::endl;
            
           // std::cout << "In triangle: " << isInTriangle(newX, newY, x1, x2, y2) << " Far enough: " << isFarEnoughFromExistingPoints(newX, newY, grid, minDistance, cellSize) << std::endl;

            if (isInTriangle(newX, newY, x1, x2, y2) &&
                isFarEnoughFromExistingPoints(newX, newY, grid, minDistance, cellSize)) {
                activePoints.push(Point(newX, newY, 0));
                points.push_back(Point(newX, newY, 0));
                grid[(int)(newX / cellSize)][(int)(newY / cellSize)] = true;
            }
        }
    }

    return points;
}


*/


double geodesicDistancePoints(Surface_mesh mesh, const Point source, const Point target){
    
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
    /*
    Print_visitor print_visitor(mesh);
    for (size_t i = 0; i < sequence_collector.sequence.size(); ++i)
      boost::apply_visitor(print_visitor, sequence_collector.sequence[i]);
    */
    double dist = 0;
    Point  prePoint = points_in_path[0];
    
    for (int i = 1; i < points_in_path.size(); ++i){
        dist = dist + sqrt((points_in_path[i].x()-prePoint.x())*(points_in_path[i].x()-prePoint.x()) + (points_in_path[i].y()-prePoint.y())*(points_in_path[i].y()-prePoint.y()) + (points_in_path[i].z()-prePoint.z())*(points_in_path[i].z()-prePoint.z()));
        prePoint = points_in_path[i];
    }
    
    return dist;
    
}

std::vector<std::vector<std::vector<bool>>> addPointToGrid(Point point, std::vector<std::vector<std::vector<bool>>> grid, double minX, double maxX, double minY, double maxY, double minZ, double maxZ, double cellSize){
    
    int px =(1 / cellSize)*point.x() - minX*(1 / cellSize);
    int py =(1 / cellSize)*point.y() - minY*(1 / cellSize);
    int pz =(1 / cellSize)*point.z() - minZ*(1 / cellSize);
    
    grid[px][py][pz] = true;
    
    
    return grid;
    
}

// Function to check if a point is within the minimum distance of existing points
bool isFarEnoughFromExistingPoints(Point point, const std::vector<std::vector<std::vector<bool>>>& grid, double minX, double maxX, double minY, double maxY, double minZ, double maxZ, double minDistance, double cellSize) {
    
    int gridX = (1 / cellSize)*point.x() - minX*(1 / cellSize);
    int gridY = (1 / cellSize)*point.y() - minY*(1 / cellSize);
    int gridZ = (1 / cellSize)*point.z() - minZ*(1 / cellSize);

    int minCheckX = fmax(0, gridX - 2);
    int maxCheckX = fmin((int)grid.size() - 1, gridX + 2);
    int minCheckY = fmax(0, gridY - 2);
    int maxCheckY = fmin((int)grid[0].size() - 1, gridY + 2);
    int minCheckZ = fmax(0, gridY - 2);
    int maxCheckZ = fmin((int)grid[0][0].size() - 1, gridZ + 2);

    for (int i = minCheckX; i <= maxCheckX; ++i) {
        for (int j = minCheckY; j <= maxCheckY; ++j) {
            for(int k = minCheckZ; k <= maxCheckZ; ++k){
                if (grid[i][j][k]) {
                    double dx = point.x() - (i * cellSize);
                    double dy = point.y() - (j * cellSize);
                    double dz = point.z() - (k * cellSize);
                    if (dx * dx + dy * dy + dz * dz < minDistance * minDistance)
                        return false;
                }
            }
        }
    }
    return true;
}


// Generate sample of random points in the annulus around
// a given point.
std::vector<Point> randomPoints(Surface_mesh mesh, Point c, int kMaxTries, double minDistance){
    
    std::vector<Point> points;
    
    
    //Begin flooding
    Face_location c_location = PMP::locate(c, mesh);
    face_descriptor fd = c_location.first;
    

    double squared_rad = 0.1 * 0.1;

    std::vector<bool> selected(num_faces(mesh), false);
    std::vector<face_descriptor> selection;
    selected[fd] = true;
    selection.push_back(fd);
    
    auto do_queue_edge = [&](halfedge_descriptor h)
    {
      halfedge_descriptor hopp=opposite(h, mesh);
      if (is_border(hopp, mesh) || selected[face(hopp, mesh)]) return false;

      K::Segment_3 edge(mesh.point(source(h,mesh)), mesh.point(target(h,mesh)));
      return K::Compare_squared_distance_3()(c, edge, squared_rad)!=CGAL::LARGER;
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

    std::cout << "center: " << c << "\n";
    for (face_descriptor f : selection)
      std::cout << f << " ";
    std::cout << "\n";
    
    Surface_mesh sub_mesh;
    
    vertex_descriptor u, v, w;
    

    halfedge_descriptor hi;
    
   
    for (face_descriptor f : selection){
        hi = mesh.halfedge(f);
     u = sub_mesh.add_vertex(mesh.point(target(hi,mesh)));
     hi = next(hi,mesh);
     v = sub_mesh.add_vertex(mesh.point(target(hi,mesh)));
     hi = next(hi,mesh);
     w = sub_mesh.add_vertex(mesh.point(target(hi,mesh)));
     sub_mesh.add_face(u,v,w);
       
    }
    
    std::cout << "Sub mesh:" << faces(sub_mesh).size() << std::endl;
    
    while(points.size() < kMaxTries){
        CGAL::Random_points_in_triangle_mesh_3<Surface_mesh>
              b(sub_mesh);
        Point d = *b;
        if(.05 < geodesicDistancePoints(mesh, d, c) && geodesicDistancePoints(mesh, d, c) < .1)
            points.push_back(d);
    }
    
    
    return points;
    
}

// Generate  Sampling Rectangle
std::vector<Point> updatedPoissonDiskSampling(Surface_mesh mesh, int sampSize, int kMaxTries, double minDistance) {
    
    std::vector<Point> points;
    std::queue<Point> activePoints;
    
    CGAL::Random_points_in_triangle_mesh_3<Surface_mesh>
          g(mesh);
    std::copy_n( g, 1, std::back_inserter(points));
    activePoints.push(points[0]);
    
    return points;
    
}

int main(int argc, char* argv[])
{
    
    const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("../eight.off");
    Surface_mesh mesh;
    if(!PMP::IO::read_polygon_mesh(filename, mesh))
    {
        std::cerr << "Invalid input." << std::endl;
        return 1;
    }
    /*
    auto vnormals = mesh.add_property_map<vertex_descriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;
    auto fnormals = mesh.add_property_map<face_descriptor, Vector>("f:normals", CGAL::NULL_VECTOR).first;
    PMP::compute_normals(mesh, vnormals, fnormals);
    std::cout << "Vertex normals :" << std::endl;
    for(vertex_descriptor vd: vertices(mesh))
        std::cout << vnormals[vd] << std::endl;
    std::cout << "Face normals :" << std::endl;
    for(face_descriptor fd: faces(mesh))
        std::cout << fnormals[fd] << std::endl;
    
    

    std::vector<Point> pointz;
    CGAL::Polygon_mesh_processing::sample_triangle_mesh(mesh, std::back_inserter(pointz), CGAL::parameters::use_grid_sampling(true).grid_spacing(.5));

    std::cout << "Sampling output :" << std::endl;
    std::cout << "******************" << std::endl;
    for (const auto& point : pointz) {
       std::cout <<  point.x() << " " <<  point.y() << " " <<  point.z()  << std::endl;
    }
     */
    std::cout << "Random points :" << std::endl;
    std::cout << "******************" << std::endl;
    std::vector<Point> points;
  
    CGAL::Random_points_in_triangle_mesh_3<Surface_mesh>
          g(mesh);
     // Get 100 random points in cdt
     std::copy_n( g, 100, std::back_inserter(points));
     // Check that we have really created 100 points.
     assert( points.size() == 100);
     // test if the generated points are close
     std::cout << points[0] << std::endl;
    
   
    
    
    std::vector<Point> pointzz = updatedPoissonDiskSampling(mesh,7,kMaxTries,kMinDistance);
    std::cout << "Function output :" << std::endl;
    std::cout << "******************" << std::endl;
    for (const auto& point : pointzz) {
       std::cout <<  point.x() << " " <<  point.y() << " " <<  point.z()  << std::endl;
    }
    
    


    /*

   
    // Triangle
    
    Point z2(100, 200, 300), z1(200, 500, 0), z3(300, 100, 600);
    
    
    std::vector<Point> pointz = poissonDiskSamplingOnTriangle(z3,z1,z2);
   
    std::cout << "Sampling output :" << std::endl;
    for (const auto& point : pointz) {
       std::cout << "(" << point.x() << ", " << point.y() << ", " << point.z() << ")" << std::endl;
    }
    
    std::cout << "Input mesh has " << faces(mesh).size() << " faces" << std::endl;
    
    
    BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
       std::cout <<  "        " << mesh.point(vd) << "\n";
     }
    
    BOOST_FOREACH(face_descriptor fd, faces(mesh)){
        std::cout <<  "        " << mesh.point(target(mesh.halfedge(fd),mesh)) << "\n";
      }
     */
    

    

    
    std::cout << "Distance Test:" << geodesicDistancePoints(mesh, points[0],points[1]) << std::endl;
    
    
    
    
    //construct a grid
    double cellSize = .1 / sqrt(2.0);
    
    
    
    double maxX = -100;
    double minX = 100;
    double maxY = -100;
    double minY = 100;
    double maxZ = -100;
    double minZ = 100;
    
    BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
        if(mesh.point(vd).x() > maxX)
            maxX = mesh.point(vd).x();
        if(mesh.point(vd).x() < minX)
            minX = mesh.point(vd).x();
        if(mesh.point(vd).y() > maxY)
            maxY = mesh.point(vd).y();
        if(mesh.point(vd).y() < minY)
            minY = mesh.point(vd).y();
        if(mesh.point(vd).z() > maxZ)
            maxZ = mesh.point(vd).z();
        if(mesh.point(vd).z() < minZ)
            minZ = mesh.point(vd).z();

     }

    int gridX = (maxX - minX) / cellSize + 1;
    int gridY = (maxY - minY) / cellSize + 1;
    int gridZ = (maxZ - minZ) / cellSize + 1;

    
    std::vector<std::vector<std::vector<bool>>> grid(gridX, std::vector<std::vector<bool> >(gridY, std::vector<bool>(gridZ)));
    
    Point c = *g;
    grid = addPointToGrid(c, grid, minX, maxX, minY, maxY, minZ, maxZ, cellSize);

  
    
    // Printing the 3d vector
    for (int i = 0; i < grid.size(); i++) {
        for (int j = 0; j < grid[i].size(); j++) {
            for (int k = 0; k < grid[i][j].size(); k++) {
                std::cout << grid[i][j][k] << " ";
            }
            std::cout << std::endl;
        }
    }
    
  
    
    std::cout << "30 Random points on sub mesh:" << std::endl;
    std::cout << "******************" << std::endl;
    std::vector<Point> sub_points = randomPoints(mesh, c, 30, .05);
    

    
    for (const auto& point : sub_points) {
        std::cout <<  point  << std::endl;
       std::cout <<  geodesicDistancePoints(mesh, point, c)  << std::endl;
    }
  
    
    return 0;
}
