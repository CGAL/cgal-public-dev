
#ifndef CGAL_POLYNOMIAL_SAGE_TYPE_GENERATOR_H
#define CGAL_POLYNOMIAL_SAGE_TYPE_GENERATOR_H

#include <CGAL/Polynomial_traits_d.h>

namespace CGAL {

template <class T, int d>
struct Polynomial_Sage_type_generator
{
private:
  typedef typename Polynomial_Sage_type_generator<T,d-1>::Type Coeff; 
public:
  typedef CGAL::Polynomial<Coeff> Type;
};

template <class T>
struct Polynomial_Sage_type_generator<T,0>{ typedef T Type; };

} //namespace CGAL

#endif // CGAL_POLYNOMIAL_GENERATOR_H
