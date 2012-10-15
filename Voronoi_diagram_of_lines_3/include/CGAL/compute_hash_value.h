#include <CGAL/basic.h>
#include<CGAL/Nef_3/Pluecker_line_3.h>

#ifndef CGAL_COMPUTE_HASH_VALUE
#define CGAL_COMPUTE_HASH_VALUE

namespace CGAL { 
// Compute an hash value for a line base on modular arithmetic 

template <class LK_3> 
int compute_hash_value(const CGAL::Line_3<LK_3>& line_3){ 
  CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);
  int old_prime = Residue::set_current_prime(67111067);
      
  typedef  Pluecker_line_3<Cartesian_tag,LK_3> PL; 
  typedef  typename LK_3::FT FT;
  typedef  Fraction_traits<FT> FTraits; 
  typedef  typename FTraits::Numerator_type NT; 
  typename FTraits::Decompose decompose;  

  PL pl(line_3);
  pl.normalize();
  CGAL::Residue r(0);
  for(int i = 0; i <6 ;i++){
    NT num,den;
    decompose(pl[i],num,den);
    r += CGAL::ipower(modular_image(num)+5,i+1) + CGAL::ipower(modular_image(den)+8,i+10);
  }
  int hash_value = r.get_value();
  if (hash_value < 0) hash_value +=  r.get_current_prime(); 
  Residue::set_current_prime(old_prime);
  return hash_value; 
}

} // namespace CGAL
#endif // CGAL_COMPUTE_HASH_VALUE
