///This is a simple c++ program for CGAL task

/*
Inplement requirement:
1. Give a convex polygon and define some particular function values at the vertices. 
2. Generate several random points inside this polygon.
3. Inplement three generic barycentric coordinate methods WP&MV&DH, compute coordinates for each inside generated points.
4. Output the interpolated function values at each inside points.
*/

/*
Notice: this simple program is generated based on the examples and tests from CGAL manuals <href:= https://doc.cgal.org/latest/Barycentric_coordinates_2/index.html>
*/

#include <iostream>
#include <math.h>
//CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Barycentric_coordinates_2/Generalized_barycentric_coordinates_2.h>

#include <CGAL/Barycentric_coordinates_2/Wachspress_2.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_2.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_2.h>

// Some(A Lot of) convenient typedefs.
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT      Scalar;
typedef Kernel::Point_2 Point;
typedef std::vector<Scalar> Scalar_vector;
typedef std::vector<Point>  Point_vector;
typedef CGAL::Creator_uniform_2<double, Point> Creator;
typedef std::back_insert_iterator<Scalar_vector> Vector_insert_iterator;

typedef CGAL::Barycentric_coordinates::Wachspress_2<Kernel> Wachspress;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Wachspress, Kernel> Wachspress_coordinates;

typedef CGAL::Barycentric_coordinates::Discrete_harmonic_2<Kernel> Discrete_harmonic;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Discrete_harmonic, Kernel> Discrete_harmonic_coordinates;

typedef CGAL::Barycentric_coordinates::Mean_value_2<Kernel> Mean_value;
typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Mean_value, Kernel> Mean_value_coordinates;

typedef boost::optional<Vector_insert_iterator> Output_type;

using std::cout; 
using std::endl; 
using std::string;
using namespace std;

int main(){
    //-1-Construct a strictly convex polygon
    const int vertice_number=12;
    Point_vector vertices(vertice_number);
    const double pi=3.14159265358;
    cout<<"strictly convex polygon:"<<endl;
    for(int i=0;i<vertice_number;i++){
        vertices[i]=Point(2*cos(i*2*pi/12.0),2*sin(i*2*pi/12.0));
        cout<<"vertice "<<i+1<<' '<<vertices[i].x()<<' '<<vertices[i].y()<<endl;
    }

    //-2-Instantiate the class with WP&DH&MV coordinates for the polygon defined above.
    Discrete_harmonic_coordinates discrete_harmonic_coordinates(vertices.begin(), vertices.end());
    Wachspress_coordinates wachspress_coordinates(vertices.begin(), vertices.end());
    Mean_value_coordinates mean_value_coordinates(vertices.begin(), vertices.end());

    //-3-Instantiate the random interior points set
    const int random_interior_points_number=100;
    Point_vector points;
    CGAL::Random_points_in_square_2<Point,Creator> point_generator(1.0);
    CGAL::cpp11::copy_n(point_generator, random_interior_points_number, std::back_inserter(points));

    //-4-Define special function values at each polygon vertices
    double defined_functions[12];
    for(int i=0;i<vertice_number;i++){
        double i_values=vertices[i].x()*vertices[i].y();
        defined_functions[i]=i_values;
    }

    //-5-Users' choice part (WP&DH&MV)
    int user_choice;
    cout<<endl<<"Please choose a particular barycentric coordinates method:"<<endl<<"1. Wachspress"<<endl<<"2. DiscreteHarmonic"<<endl<<"3. MeanValue"<<endl;
    cout<<"Press 1 or 2 or 3 to choose one methods"<<endl<<endl;
    cin>>user_choice;

    //-6-Compute coordinates and target function values at each interior points based on the chosen method
    switch(user_choice){
        case 1:
        cout << endl << "Computed Wachspress coordinates: " << endl << endl;
        for(int i = 0; i < random_interior_points_number; ++i) {
            // Compute coordinates.
            Scalar_vector coordinates;
            coordinates.reserve(vertice_number);
            wachspress_coordinates(points[i], std::back_inserter(coordinates));
            //Compute target function values
            double target_values=0.0;
            for(int j=0;j<vertice_number;j++){
                target_values+=coordinates[j]*defined_functions[j];
            }
            //Output values
            cout<<"Point "<<i+1<<": "<<target_values<<endl; 
        }
        break;
        
        case 2:
        cout << endl << "Computed DiscreteHarmonic coordinates: " << endl << endl;
        for(int i = 0; i < random_interior_points_number; ++i) {
            // Compute coordinates.
            Scalar_vector coordinates;
            coordinates.reserve(vertice_number);
            discrete_harmonic_coordinates(points[i], std::back_inserter(coordinates));
            //Compute target function values
            double target_values=0.0;
            for(int j=0;j<vertice_number;j++){
                target_values+=coordinates[j]*defined_functions[j];
            }
            //Output values
            cout<<"Point "<<i+1<<": "<<target_values<<endl; 
        }
        break;
        
        case 3:
        cout << endl << "Computed MeanValue coordinates: " << endl << endl;
        for(int i = 0; i < random_interior_points_number; ++i) {
            // Compute coordinates.
            Scalar_vector coordinates;
            coordinates.reserve(vertice_number);
            mean_value_coordinates(points[i], std::back_inserter(coordinates));
            //Compute target function values
            double target_values=0.0;
            for(int j=0;j<vertice_number;j++){
                target_values+=coordinates[j]*defined_functions[j];
            }
            //Output values
            cout<<"Point "<<i+1<<": "<<target_values<<endl; 
        }
        break;
        
        default:
        break;
    }

    return 0;
}

