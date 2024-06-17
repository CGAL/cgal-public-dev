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

//Function to generate a random point in a triangle
Point randomInTriangle(double x1, double x2, double y2){
    
    double u = (double)rand() / RAND_MAX;
    double v = (double)rand() / RAND_MAX;
    
    // Make sure u + v <= 1
        if (u + v > 1.0) {
            u = 1.0 - u;
            v = 1.0 - v;
        }
    
    double w = 1.0 - u - v;
    Point randomPoint(u * x1 + v * x2,u * 0 + v * y2,0);
    
    return randomPoint;
}

// Function to check if a point is in a rectangle
bool isInRectangle(double x, double y, double rectWidth, double rectHeight) {
    return x >= 0 && x <= rectWidth && y >= 0 && y <= rectHeight;
}

// Function to check if a point is in a triangle
bool isInTriangle(double x, double y, double x1, double x2, double y2) {
    return y >= 0 && y <= (y2/x2)*x && y <= (y2/(x2-x1))*(x-x1);
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

//Normalize a point
Point normalize(Point p){
    double length = sqrt(p.x()*p.x()+p.y()*p.y()+p.z()*p.z());
    Point n(p.x()/length, p.y()/length, p.z()/length);
    return n;
}
//Normalized normal vector of two points
Point crossProduct(Point p1, Point p2){
    
    Point normal(p1.y()*p2.z()-p1.z()*p2.y(),-1*(p1.x()*p2.z()-p1.z()*p2.x()), p1.x()*p2.y()-p1.y()*p2.x());
    
    return normal;
}



//A function to Poisson sample on a triangle
std::vector<Point> poissonDiskSamplingOnTriangle(Point q1, Point q2, Point q3){
    
    //Make p1 p2 the longest side of the triangle
    double d12 = (q2.x()-q1.x())*(q2.x()-q1.x())+(q2.y()-q1.y())*(q2.y()-q1.y())+(q2.z()-q1.z())*(q2.z()-q1.z());
    double d13 = (q3.x()-q1.x())*(q3.x()-q1.x())+(q3.y()-q1.y())*(q3.y()-q1.y())+(q3.z()-q1.z())*(q3.z()-q1.z());
    double d23 = (q3.x()-q2.x())*(q3.x()-q2.x())+(q3.y()-q2.y())*(q3.y()-q2.y())+(q3.z()-q2.z())*(q3.z()-q2.z());
    
    double max_edge_squared = d12;
    Point p1 = q1;
    Point p2 = q2;
    Point p3 = q3;
    
    if(d13>max_edge_squared){
        max_edge_squared = d13;
        p1 = q1;
        p2 = q3;
        p3 = q2;
    }
    if(d23>max_edge_squared){
        max_edge_squared = d23;
        p1 = q2;
        p2 = q3;
        p3 = q1;
    }
    
    //Find rotation matrix so that triangle has an edge on x-axis
    Point p4(p2.x()-p1.x(), p2.y()-p1.y(),p2.z()-p1.z());
    Point U(normalize(p4));
    
 
    Point p5(p3.x()-p1.x(), p3.y()-p1.y(),p3.z()-p1.z());
    Point W(normalize(crossProduct(U,p5)));
    
    Point V(normalize(crossProduct(U,W)));
    
    std::cout << "(" << U.x() << ", " << U.y() << ", " << U.z() << ")" << std::endl;
    std::cout << "(" << V.x() << ", " << V.y() << ", " << V.z() << ")" << std::endl;
    std::cout << "(" << W.x() << ", " << W.y() << ", " << W.z() << ")" << std::endl;
    
    //The matrix [U,V,X]^T maps points in U V space to standard basis
    //The matrix [U,V,W] maps the standard basis to U V space
    //The shifted triangle has points (0,0,0), (x1, 0, 0), and (x2, y2, 0)
    double x1 = U.x()*p4.x()+U.y()*p4.y()+U.z()*p4.z();
    double x2 = U.x()*p5.x()+U.y()*p5.y()+U.z()*p5.z();
    double y2 = V.x()*p5.x()+V.y()*p5.y()+V.z()*p5.z();
    

    //Note y2 could still be negative
    std::cout << "(" << x1 << "," << x2 << ", " << y2 <<  ")" << std::endl;
    
    //Call generate Poisson Disk Sample
    std::vector<Point> pointsInitial = generatePoissonDiskSampling(x1, x2, abs(y2), kMinDistance);
    
    //Put points back in original triangle
    //Flip y-coordinate if y2 is negative
    std::vector<Point> pointsUpdate;
    if(y2 < 0){
        for (const auto& point : pointsInitial) {
            Point update(U.x()*point.x()+V.x()*point.y()+W.x()*point.z()+p1.x(), -1*(U.y()*point.x()+V.y()*point.y()+W.y()*point.z())+p1.y(), U.z()*point.x()+V.z()*point.y()+W.z()*point.z()+p1.z());
            pointsUpdate.push_back(update);
        }
    }else{
        for (const auto& point : pointsInitial) {
            Point update(U.x()*point.x()+V.x()*point.y()+W.x()*point.z()+p1.x(), U.y()*point.x()+V.y()*point.y()+W.y()*point.z()+p1.y(), U.z()*point.x()+V.z()*point.y()+W.z()*point.z()+p1.z());
            pointsUpdate.push_back(update);
        }
    }
    
    return pointsUpdate;
}




int main(int argc, char* argv[])
{
    /*
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

    
    
    
    // Rectangle
    double width = 800;
    double height = 600;
    
    std::vector<Point> points = generatePoissonDiskSampling(width, height, kMinDistance);
    
    
    for (const auto& point : points) {
        std::cout << "(" << point.x() << ", " << point.y() << ", " << point.z() << ")" << std::endl;
    }
    */
    // Triangle
    
    Point z2(100, 200, 300), z1(200, 500, 0), z3(300, 100, 600);
    
    
    std::vector<Point> pointz = poissonDiskSamplingOnTriangle(z3,z1,z2);
   
    std::cout << "Sampling output :" << std::endl;
    for (const auto& point : pointz) {
       std::cout << "(" << point.x() << ", " << point.y() << ", " << point.z() << ")" << std::endl;
    }

    return 0;
}
