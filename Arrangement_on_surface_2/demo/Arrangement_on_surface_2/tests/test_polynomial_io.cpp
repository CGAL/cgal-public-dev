#include "PolynomialParser.h"

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

Polynomial_2 Load( const std::string& fn )
{
  Polynomial_2 res;
  std::ifstream ifs( fn.c_str( ) );
  ifs >> res;
  ifs.close( );
  return res;
}

////////////////////////////////////////////////////////////////////////////
//  Main program
////////////////////////////////////////////////////////////////////////////
int
main()
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
