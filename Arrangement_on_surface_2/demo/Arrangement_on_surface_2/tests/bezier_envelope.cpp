#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
#error "CGAL needs CORE."
#endif

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_default_dcel.h>

#include <CGAL/envelope_2.h>
#include <CGAL/Envelope_diagram_1.h>

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
  Bezier_traits_2;
typedef Bezier_traits_2::Curve_2 Bezier_curve_2;
typedef Bezier_traits_2::X_monotone_curve_2 Bezier_x_monotone_curve_2;
typedef CGAL::Arrangement_with_history_2<Bezier_traits_2>
  Bezier_arrangement_2;
typedef CGAL::Envelope_diagram_1< Bezier_traits_2 > Diagram_1;

int LoadArr( const std::string& fn,
  Bezier_arrangement_2* arr )
{
  // Get the name of the input file from the command line, or use the default
  // Bezier.dat file if no command-line parameters are given.

  // Open the input file.
  std::ifstream in_file(fn.c_str( ));

  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << fn << std::endl;
    return 1;
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
  insert (*arr, curves.begin(), curves.end());

  // Print the arrangement size.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr->number_of_vertices()
            << ",  E = " << arr->number_of_edges()
            << ",  F = " << arr->number_of_faces() << std::endl;
  return 0;
}

void BuildEnvelope( Bezier_arrangement_2& arr, bool upper,
  Diagram_1* diagram )
{
  std::list< Bezier_x_monotone_curve_2 > curves;
  Bezier_arrangement_2::Edge_iterator eit;
  for (eit = arr.edges_begin( ); eit != arr.edges_end( ); ++eit)
  {
    curves.push_back( eit->curve( ) );
  }

  if ( ! upper )
  {
    CGAL::lower_envelope_x_monotone_2(curves.begin(), curves.end(), *diagram);
  }
  else
  {
    CGAL::upper_envelope_x_monotone_2(curves.begin(), curves.end(), *diagram);
  }
}

int main( int argc, char *argv[] )
{
  std::string input = "Bezier.dat";
  if ( argc >= 2 )
  {
    input = argv[1];
  }
  else
  {
    std::cout << "Looking in current directory for Bezier.dat...\n";
  }

  Bezier_arrangement_2 arr;
  if ( LoadArr( input, &arr ) )
  {
    std::cout << "Error loading data from " << input << "\n";
    return 1;
  }

  Diagram_1 upper_envelope;
  Diagram_1 lower_envelope;
  std::cout << "Constructing lower envelope...\n";
  BuildEnvelope( arr, false, &lower_envelope );
  std::cout << "Constructing upper envelope...\n";
  BuildEnvelope( arr, true, &upper_envelope );

  return 0;
}
