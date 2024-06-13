#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
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



int kMaxTries = 30; // Number of attempts to find a point
double kMinDistance = 20.0; // Minimum distance between points

// Function to check if a point is in a rectangle
bool isInRectangle(double x, double y, double rectWidth, double rectHeight) {
    return x >= 0 && x <= rectWidth && y >= 0 && y <= rectHeight;
}


// Function to check if a point is within the minimum distance of existing points
bool isFarEnoughFromExistingPoints(double x, double y, const std::vector<std::vector<bool>>& grid, double minDistance, double cellSize) {
    int gridX = x / cellSize;
    int gridY = y / cellSize;

    int minX = fmax(0, gridX - 2);
    int maxX = fmin((int)grid.size() - 1, gridX + 2);
    int minY = fmax(0, gridY - 2);
    int maxY = fmin((int)grid[0].size() - 1, gridY + 2);

    for (int i = minX; i <= maxX; ++i) {
        for (int j = minY; j <= maxY; ++j) {
            if (grid[i][j]) {
                double dx = x - (i * cellSize);
                double dy = y - (j * cellSize);
                if (dx * dx + dy * dy < minDistance * minDistance)
                    return false;
            }
        }
    }
    return true;
}


// Generate  Sampling
std::vector<Point> generatePoissonDiskSampling(double width, double height, double minDistance) {
    srand(time(NULL));
    double cellSize = minDistance / sqrt(2.0);

    int gridWidth = (int)ceil(width / cellSize);
    int gridHeight = (int)ceil(height / cellSize);

    std::vector<std::vector<bool>> grid(gridWidth, std::vector<bool>(gridHeight, false));
    std::vector<Point> points;
    std::queue<Point> activePoints;

    // Generate first point
    double startX = width * (double)rand() / RAND_MAX;
    double startY = height * (double)rand() / RAND_MAX;
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

            if (isInRectangle(newX, newY, width, height) &&
                isFarEnoughFromExistingPoints(newX, newY, grid, minDistance, cellSize)) {
                activePoints.push(Point(newX, newY, 0));
                points.push_back(Point(newX, newY, 0));
                grid[(int)(newX / cellSize)][(int)(newY / cellSize)] = true;
            }
        }
    }

    return points;
}

//Transform triangle in three space to have a vertex at the origin and another
//vertex on the x-axis
std::vector<Point> affineTransformTriangle(Point p1, Point p2, Point p3){
    
    double matrix[3][3];
    matrix[0][0]=0;
    matrix[0][1]=0;
    matrix[0][2]=0;
    
    matrix[1][0]=p2.x()-p1.x();
    matrix[1][1]=p2.y()-p1.y();
    matrix[1][2]=p2.z()-p1.z();
    
    matrix[2][0]=p3.x()-p1.x();
    matrix[2][1]=p3.y()-p1.y();
    matrix[2][2]=p3.z()-p1.z();

    double normp2 = sqrt(matrix[1][0]*matrix[1][0]+matrix[1][1]*matrix[1][1]+matrix[1][2]*matrix[1][2]);
    double sq_dist_p1p3 = matrix[2][0]*matrix[2][0]+matrix[2][1]*matrix[2][1]+matrix[2][2]*matrix[2][2];
    double sq_dist_p2p3 = pow(matrix[2][0]-matrix[1][0],2)+pow(matrix[2][1]-matrix[1][1],2)+pow(matrix[2][2]-matrix[1][2],2);
    
    
    Point p4(0, 0, 0), p5(normp2, 0, 0), p6((sq_dist_p2p3-sq_dist_p1p3-normp2*normp2)/(-2*normp2), sqrt(sq_dist_p1p3-matrix[2][0]*matrix[2][0]), 0);
    
    std::vector<Point> points;
    points.push_back(p4);
    points.push_back(p5);
    points.push_back(p6);
    
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
    auto vnormals = mesh.add_property_map<vertex_descriptor, Vector>("v:normals", CGAL::NULL_VECTOR).first;
    auto fnormals = mesh.add_property_map<face_descriptor, Vector>("f:normals", CGAL::NULL_VECTOR).first;
    PMP::compute_normals(mesh, vnormals, fnormals);
    std::cout << "Vertex normals :" << std::endl;
    for(vertex_descriptor vd: vertices(mesh))
        std::cout << vnormals[vd] << std::endl;
    std::cout << "Face normals :" << std::endl;
    for(face_descriptor fd: faces(mesh))
        std::cout << fnormals[fd] << std::endl;
    
    
    // New things

    std::cout << "Sampling output :" << std::endl;
    
    
    // Rectangle
    double width = 800;
    double height = 600;
    
    std::vector<Point> points = generatePoissonDiskSampling(width, height, kMinDistance);
    
    
    for (const auto& point : points) {
        std::cout << "(" << point.x() << ", " << point.y() << ", " << point.z() << ")" << std::endl;
    }
    
    Point p1(1, 3, 4), p2(4, 1, 3), p3(5, 2, 1);
    
    for (const auto& point : affineTransformTriangle(p1, p2, p3)) {
        std::cout << "(" << point.x() << ", " << point.y() << ", " << point.z() << ")" << std::endl;
    }
    
   

    return 0;
}
