#include "PolynomialParser.h"
#include "AlgebraicDemoTraits.h"
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QApplication>
#include "Utils.h"

typedef AlgebraicDemoTraits::ArrTraitsType ArrTraitsType;
typedef ArrTraitsType::Polynomial_2 Polynomial_2;
typedef PolynomialParser< Polynomial_2 > ParserType;


Polynomial_2 Parabola( )
{
  std::string str = "y - x^2";
  Polynomial_2 parabola;
  ParserType::Parse( str, &parabola );
  return parabola;
}

Polynomial_2 Degree6Curve( )
{
  std::string str = "x^6 + y^6 - x^3y^3 - 12";
  Polynomial_2 res;
  ParserType::Parse( str, &res );
  return res;
}

template < typename OutputIterator >
int Load( const std::string& fn, OutputIterator oit )
{
  int count = 0;
  std::ifstream ifs( fn.c_str( ) );
  std::string str;
  while ( std::getline( ifs, str ) )
  {
    Polynomial_2 polynomial;
    bool ok = ParserType::Parse( str, &polynomial );
    if ( ok )
    {
      *oit = polynomial;
      ++oit;
      ++count;
    }
  }
  ifs.close( );
  return count;
}

int TestParser( )
{
  CGAL::set_pretty_mode( std::cout );
  std::cout << "/////////////////////////////////////////////////////////\n\n";
  std::cout << "\t\tA polynomial parser for Spirit...\n\n";
  std::cout << "/////////////////////////////////////////////////////////\n\n";

  std::cout
    << "Give me a polynomial with terms of the form :"
    << "\"A*x^Ny^M\", where A, N, and M are integers\n";
  std::cout << "Type [q or Q] to quit\n\n";

  std::string str;
  while (getline(std::cin, str))
  {
    if (str.empty() || str[0] == 'q' || str[0] == 'Q')
      break;

    Polynomial_2 po;
    bool ok = ParserType::Parse( str, &po );

    if ( ok )
    {
      std::cout << "-------------------------\n";
      std::cout << "Parsing succeeded\n";
      std::cout << po << "\n";
      std::cout << "\n-------------------------\n";
    }
    else
    {
      std::cout << "-------------------------\n";
      std::cout << "Parsing failed\n";
      std::cout << "-------------------------\n";
    }
  }

  std::cout << "Bye... :-) \n\n";
  return 0;
}

int main( int argc, char *argv[] )
{
  typedef AlgebraicDemoTraits::ArrTraitsType TraitsType;
  typedef AlgebraicDemoTraits::ArrangementType ArrangementType;
  typedef TraitsType::Curve_2 Curve_2;
  typedef ArrangementType::Curve_iterator Curve_iterator;
  typedef TraitsType::X_monotone_curve_2 X_monotone_curve_2;
  TraitsType traits;

  CGAL::set_pretty_mode( std::cout );
  if ( argc < 2 )
  {
    std::cout << "Usage: " << argv[0] << " polynomials-files\n";
    return 0;
  }

  ArrangementType arr;

  std::vector< Polynomial_2 > polynomials;
  int count = Load( argv[1], std::back_inserter( polynomials ) );
  std::cout << count << " polynomials loaded\n";

  for ( int i = 0; i < polynomials.size( ); ++i )
  {
    std::cout << polynomials[i] << "\n";
  }

  TraitsType::Construct_curve_2 construct_curve =
    traits.construct_curve_2_object( );
  std::vector< Curve_2 > curves;
  for ( int i = 0; i < polynomials.size( ); ++i )
  {
    Curve_2 cv = construct_curve( polynomials[i] );
    curves.push_back( cv );
  }

  CGAL::insert( arr, curves.begin( ), curves.end( ) );

  std::cout << "The arrangement size:" << std::endl
    << " V = " << arr.number_of_vertices()
    << ", E = " << arr.number_of_edges()
    << ", F = " << arr.number_of_faces() << std::endl;

  TraitsType::Make_x_monotone_2 make_x_monotone_2 =
    traits.make_x_monotone_2_object( );

  for ( Curve_iterator it = arr.curves_begin( );
        it != arr.curves_end( );
        ++it )
  {
    std::vector< CGAL::Object > cvs;
    make_x_monotone_2( *it, std::back_inserter( cvs ) );
    for ( unsigned int i = 0 ; i < cvs.size( ); ++i )
    {
      X_monotone_curve_2 cv;
      CGAL::assign( cv, cvs[ i ] );
      CGAL::Bbox_2 bb = cv.bbox( );
    }
  }

  return 0;
}
