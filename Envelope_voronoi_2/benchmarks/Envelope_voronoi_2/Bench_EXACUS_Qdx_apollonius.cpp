#include <CGAL/EXACUS_apollonius_traits_2.h>
#include "bench_classes.h"

// EXACUS apollonius using quadrics
void Bench_EXACUS_Qdx_apollonius::convert_sites_to_quadrics()
{
  typedef AT::Poly_rat3 Poly_rat3;

  Base_sites::iterator it;
  for (it = _sites.begin(); it != _sites.end(); ++it)
  {
    Rational x = it->center().first;
    Rational y = it->center().second;
    Rational r = -it->r();
    
    // construct polynomial (x-d[0])^2
    Poly_rat3 pr = QdX::get_polynomial_from_coefficients(
      Rational(1), Rational(0), Rational(0),
      Rational(1), Rational(0),
      Rational(-1),
      Rational(-2*x), Rational(-2*y), Rational(2*r),
      Rational((x*x) + (y*y) - (r*r))
      );
    
    typedef CGAL::Fraction_traits< Poly_rat3 > FT;
    FT::Numerator_type pi;
    FT::Denominator_type dummy;
    FT::Decompose decompose;
    decompose(pr,pi,dummy);

    P_quadric_3 q = P_quadric_3(pi);
    
    std::cerr << "Polinomial for Quadric: " << pr << std::endl;
    
    _quadric_sites.push_back(q);
  }
  
  
  
}

int Bench_EXACUS_Qdx_apollonius::init()
{
  Base::init();
  ::CGAL::set_mode(std::cout, ::CGAL::IO::PRETTY);
  ::CGAL::set_mode(std::cerr, ::CGAL::IO::PRETTY);

  convert_sites_to_quadrics();

  return 0;
}

void Bench_EXACUS_Qdx_apollonius::op()
{
  EXACUS_Qdx_apollonius_diagram diagram;
  CGAL::lower_envelope_3(_quadric_sites.begin(), _quadric_sites.end(), diagram);
  
  if (m_verbose_level > 0)
  {
    std::cout << "# of vertices: " << diagram.number_of_vertices() << std::endl;
    std::cout << "# of halfedges: " << diagram.number_of_halfedges() << std::endl;
    std::cout << "# of faces: " << diagram.number_of_faces() << std::endl;
    
    if (m_verbose_level > 1) 
    {
      EXACUS_Qdx_apollonius_diagram::Vertex_iterator it = 
        diagram.vertices_begin();
      for (; it != diagram.vertices_end(); ++it)
      {
        std::cout << it->point() << std::endl;
      }
    }
  }
}

