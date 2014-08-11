#include <QObject>
#include <CGAL/Simple_cartesian.h>
#include "ArrangementDemoTraits.h"
#include "Utils.h"

typedef BezierDemoTraits::ArrangementType ArrType;
typedef ArrType::Curve_2 Curve_2;
typedef ArrType::X_monotone_curve_2 X_monotone_curve_2;
typedef CGAL::Simple_cartesian<double> App_kernel;
typedef App_kernel::Point_2 App_point_2;
typedef App_kernel::Segment_2 App_segment_2;

template < class OutputIterator >
void get_control_pts( const X_monotone_curve_2& curve,
  OutputIterator oit )
{
    Curve_2 supportingCurve = curve.supporting_curve();
    for (int k = 0; k < supportingCurve.number_of_control_points(); k++)
    {
      const Rat_point_2& pt = supportingCurve.control_point(k);
      *oit = App_point_2 (CGAL::to_double (pt.x()),
                                     CGAL::to_double (pt.y()));
      ++oit;
    }
}

double estimate_alpha( const X_monotone_curve_2& curve )
{
  std::vector< std::pair<double, double> > samples;
  std::pair<double, double> range = curve.parameter_range( );
  Curve_2 supportingCurve = curve.supporting_curve();
  int P = 2 * (supportingCurve.number_of_control_points() + 1);
  supportingCurve.sample(range.first, range.second,
    2*P + 1, // number of samples
    std::back_inserter( samples ) );

  std::vector< App_point_2 > query_pts;
  std::vector< App_segment_2 > segments;
  double max_error = 0.0;
  for (int i = 0; i < P; ++i)
  {
    App_point_2 p1( samples[2*i].first, samples[2*i].second );
    App_point_2 p2( samples[2*i+2].first, samples[2*i+2].second );
    App_segment_2 seg( p1, p2 );
    segments.push_back( seg );
  }

  for (int i = 0; i < P; ++i)
  {
    App_point_2 query_pt( samples[2*i+1].first, samples[2*i+1].second );
    double error = CGAL::squared_distance( query_pt, segments[i] );
    if ( error > max_error )
      max_error = error;
  }

  return sqrt(max_error);
}

double test( X_monotone_curve_2& curve )
{
  std::vector< App_point_2 > ctrl_pts;
  get_control_pts( curve, std::back_inserter( ctrl_pts ) );
  Curve_2 supportingCurve = curve.supporting_curve();
  double P = 2 * (supportingCurve.number_of_control_points() + 1);
  // TODO: Calculate this number appropriately
  int N = P;

  double alpha = estimate_alpha( curve );
  const double Epsilon = 0.01;
  const int Factor = P / 2 * sqrt( alpha / Epsilon );
  N = ( Factor > N )? Factor : N;
  std::cout << "alpha: " << alpha << "\n";
  std::cout << "N: " << N << "\n";

  // TODO: calculate and cache alpha for each curve seen
  // int N = P / 2 * sqrt(alpha / epsilon)
  // alpha is the max chordal deviation for sampling width 1/P
  // epsilon is the desired error bound
  // TODO: Scale the number of samples by the norm of the view matrix
  std::vector< App_point_2 > samples;

  return 0.0;
}

int main( int argc, char *argv[] )
{
  ArrType arr;
  LoadArrFromFile< BezierDemoTraits > loadArr;
  loadArr( argv[1], &arr );
  test( arr.edges_begin( )->curve( ) );

  return 0;
}
