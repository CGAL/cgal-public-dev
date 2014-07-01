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

/**
Draw a slash through each curve's bounding box, and save any intersection points.
\param[out] oit has a value type of std::pair< double, double >
*/
template < class OutputIterator >
void IntersectCurvePairs( Arrangement_2& arr, OutputIterator oit )
{
  typedef Arrangement_2::Curve_2 Curve_2;
  typedef Arrangement_2::X_monotone_curve_2 X_monotone_curve_2;
  typedef Arrangement_2::Edge_iterator Edge_iterator;
  typedef Arrangement_2::Geometry_traits_2 Geometry_traits_2;
  typedef Geometry_traits_2::Multiplicity Multiplicity;
  typedef Geometry_traits_2::Point_2 Point_2;
  typedef Geometry_traits_2::Intersect_2 Intersect_2;
  typedef Geometry_traits_2::Make_x_monotone_2 Make_x_monotone_2;
  typedef std::pair< Point_2, Multiplicity > IntersectionResult;

  Geometry_traits_2 traits;
  Intersect_2 intersect_curves = traits.intersect_2_object( );
  Make_x_monotone_2 make_x_monotone = traits.make_x_monotone_2_object( );

  for ( Edge_iterator it = arr.edges_begin( );
        it != arr.edges_end( ); ++it )
  {
    X_monotone_curve_2 curve = it->curve( );
    Curve_2 supporting_curve = curve.supporting_curve( );
    CGAL::Bbox_2 bbox = supporting_curve.bbox( );
    Rat_point_2 p1( bbox.xmin(), bbox.ymin() );
    Rat_point_2 p2( bbox.xmax(), bbox.ymax() );
    std::vector< Rat_point_2 > pts;
    pts.push_back( p1 );
    pts.push_back( p2 );
    Curve_2 probe_line( pts.begin( ), pts.end( ) );
    X_monotone_curve_2 x_probe_line;

    // convert probe_line to x-monotone curve
    std::vector< CGAL::Object > vo;
    make_x_monotone( probe_line, std::back_inserter( vo ) );
    assert( vo.size( ) == 1 );

    // intersect the probe line with the curve
    if ( ! CGAL::assign( x_probe_line, vo[0] ) )
    {
      continue;
    }

    CGAL::Object o;
    CGAL::Oneset_iterator< CGAL::Object > oi( o );
    intersect_curves( curve, x_probe_line, oi );

    // inspect the intersection result
    IntersectionResult pair;
    if ( CGAL::assign( pair, o ) )
    {
      Point_2 pt = pair.first;
      std::pair< double, double > approx_pt = pt.approximate();
      std::cout << "Intersection at: ("
        << approx_pt.first << ", "
        << approx_pt.second << ")\n";
      *oit = approx_pt;
      ++oit;
    }
    else
    {
      std::cout << "No intersection.\n";
    }
  }
}

int main( int argc, char* argv[] )
{
  QApplication app( argc, argv );

  // Get the name of the input file from the command line, or use the default
  // Bezier.dat file if no command-line parameters are given.
  std::string filename = (argc > 1) ? argv[1] : "Bezier.dat";

  Arrangement_2* arr = LoadArr( filename );
  std::vector< std::pair< double, double > > points;
  //IntersectCurvePairs( *arr, std::back_inserter( points ) );

  BezierExampleWindow window;
  window.setArrangement( arr );
  window.show( );
  window.drawXYPairs( points.begin(), points.end() );

  return app.exec( );
}
