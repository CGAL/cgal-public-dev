#include "BezierExampleTypes.h"
#include "BezierExampleWindow.h"
#include <QApplication>

Arrangement_2* LoadArr( const std::string& filename )
{
  // Get the name of the input file from the command line, or use the default
  // Bezier.dat file if no command-line parameters are given.

  // Open the input file.
  std::ifstream   in_file (filename.c_str( ));

  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
    return NULL;
  }

  // Read the curves from the input file.
  unsigned int               n_curves;
  std::list<Bezier_curve_2>  curves;
  Bezier_curve_2             B;
  unsigned int               k;

  in_file >> n_curves;
  for (k = 0; k < n_curves; k++) {
    // Read the current curve (specified by its control points).
    in_file >> B;
    curves.push_back (B);

    std::cout << "B = {" << B << "}" << std::endl;
  }

  // Construct the arrangement.
  Arrangement_2* arr = new Arrangement_2;
  insert (*arr, curves.begin(), curves.end());

  // Print the arrangement size.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr->number_of_vertices()
            << ",  E = " << arr->number_of_edges()
            << ",  F = " << arr->number_of_faces() << std::endl;

  return arr;
}

int main( int argc, char* argv[] )
{
  QApplication app( argc, argv );

  // Get the name of the input file from the command line, or use the default
  // Bezier.dat file if no command-line parameters are given.
  std::string filename = (argc > 1) ? argv[1] : "Bezier.dat";

  Arrangement_2* arr = LoadArr( filename );

  BezierExampleWindow window;
  window.setArrangement( arr );
  window.show( );

  return app.exec( );
}
