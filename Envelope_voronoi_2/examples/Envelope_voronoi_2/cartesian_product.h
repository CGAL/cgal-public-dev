#ifndef CGAL_CARTESIAN_PRODUCT_H
#define CGAL_CARTESIAN_PRODUCT_H

#include <CGAL/basic.h>

namespace CGAL {

template<class I, class O>
O cartesian_product(I b1, I e1, I b2, I e2, O o)
{
  for (I i = b1; i != e1; ++i)
    for (I j = b2; j != e2; ++j)
      *o++ = std::make_pair(*i, *j);
  
  return o;
}

template<class I, class O>
O strict_cartesian_product(I b1, I e1, I b2, I e2, O o)
{
  for (I i = b1; i != e1; ++i)
    for (I j = b2; j != e2; ++j)
      if (i != j)
        *o++ = std::make_pair(*i, *j);
  
  return o;
}

template<class I, class O>
O cartesian_square_without_order(I b, I e, O o)
{
  for (I i = b; i != e; ++i)
    for (I j = i; j != e; ++j)
      if (i != j)
        *o++ = std::make_pair(*i, *j);
  
  return o;
}


} //namespace CGAL

#endif // CGAL_CARTESIAN_PRODUCT_H
