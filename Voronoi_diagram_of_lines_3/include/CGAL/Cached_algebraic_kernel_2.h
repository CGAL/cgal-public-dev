// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//                 Sebastian Limbach <slimbach@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_CHACHED_ALGEBRAIC_KERNEL_2_H
#define CGAL_CHACHED_ALGEBRAIC_KERNEL_2_H

#include <CGAL/basic.h>
#include <CGAL/Cached_algebraic_kernel_1.h>
#include <map>

namespace CGAL {

// sofar this is a pure wrapper 
template< class AlgebraicKernel_2 > 
class Cached_algebraic_kernel_2
  :public Cached_algebraic_kernel_1<AlgebraicKernel_2>{
public:
  typedef Cached_algebraic_kernel_1<AlgebraicKernel_2>  Base; 
  typedef Cached_algebraic_kernel_1<AlgebraicKernel_2>  Algebraic_kernel_d_1;  

// this would be better, but is not possible due to design of AK_2 
//   :public AlgebraicKernel_2::Algebraic_kernel_d_1{       
// public:
//   typedef typename AlgebraicKernel_2::Algebraic_kernel_d_1 Algebraic_kernel_d_1;
//   typedef typename AlgebraicKernel_2::Algebraic_kernel_d_1 Base;
  
  typedef AlgebraicKernel_2 AK_2;
public:
  typedef AlgebraicKernel_2 Algebraic_kernel_d_2;
  typedef Cached_algebraic_kernel_2<Algebraic_kernel_d_2> Self;

  typedef typename Algebraic_kernel_d_2::Polynomial_2     Polynomial_2; 
  typedef typename Algebraic_kernel_d_2::Algebraic_real_2 Algebraic_real_2;
  typedef typename Algebraic_kernel_d_2::Boundary         Boundary;
  typedef typename Algebraic_kernel_d_2::Coefficient      Coefficient;   
  
#define CGAL_ALGEBRAIC_KERNEL_2_PRED(Y,Z)       \
  typedef typename Algebraic_kernel_d_2::Y Y;     \
    Y Z() const { return algebraic_kernel_2().Z();}

  // required by AlgebraicKernel_2
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Boundary_between_x_2,boundary_between_x_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Boundary_between_y_2,boundary_between_y_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Compare_x_2,compare_x_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Compare_xy_2,compare_xy_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Compare_y_2,compare_y_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Get_x_2,get_x_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Get_y_2,get_y_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Is_coprime_2,is_coprime_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Is_square_free_2,is_square_free_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Make_coprime_2,make_coprime_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Make_square_free_2,make_square_free_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Sign_at_2,sign_at_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Solve_2,solve_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(X_critical_points_2,x_critical_points_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Y_critical_points_2,y_critical_points_2_object);
  // chached, see below 
//  CGAL_ALGEBRAIC_KERNEL_2_PRED(Square_free_factorize_2,square_free_factorize_2_object);
  // deprecated 
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Upper_boundary_x_2,upper_boundary_x_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Lower_boundary_x_2,lower_boundary_x_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Upper_boundary_y_2,upper_boundary_y_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Lower_boundary_y_2,lower_boundary_y_2_object);
  
  // required by " with analyis " 
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Construct_curve_2,construct_curve_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Construct_curve_pair_2,construct_curve_pair_2_object);
  CGAL_ALGEBRAIC_KERNEL_2_PRED(Decompose_2,decompose_2_object);
#undef CGAL_ALGEBRAIC_KERNEL_2_PRED  

  // implicit required by code 
  typedef typename AlgebraicKernel_2::Coordinate_1  Coordinate_1; 
  typedef typename AlgebraicKernel_2::Coordinate_2 Coordinate_2; 
  typedef typename AlgebraicKernel_2::Curve_analysis_2 Curve_analysis_2;
  typedef typename AlgebraicKernel_2::Curve_pair_analysis_2 Curve_pair_analysis_2; 

private: 
  Algebraic_kernel_d_2 m_algebraic_kernel_2;
public:
  const Algebraic_kernel_d_2& algebraic_kernel_2() const {
    return m_algebraic_kernel_2; 
  };
  
public:
  Cached_algebraic_kernel_2(){};

  //Cached_algebraic_kernel_2(const Algebraic_kernel_d_2& ak_2)
  //  :Base(ak_2), m_algebraic_kernel_2(ak_2) {};

  virtual ~Cached_algebraic_kernel_2(){}


private:
  typedef boost::optional<std::list<std::pair<Polynomial_2,int> > > Sqff_cache_data;
  typedef std::map<Polynomial_2, Sqff_cache_data > Sqff_cache;
  mutable Sqff_cache m_sqff_cache;
public:
  struct Square_free_factorize_2{
    typedef typename AK_2::Square_free_factorize_2::first_argument_type first_argument_type;
    const Cached_algebraic_kernel_2* m_cak_2;
    Square_free_factorize_2(const Cached_algebraic_kernel_2 * cak_2):m_cak_2(cak_2){}
    template <typename OutputIterator>
    OutputIterator operator()(first_argument_type p, OutputIterator oit){
      Sqff_cache_data& sqff = m_cak_2->m_sqff_cache[p];
      if(!sqff.is_initialized()){
        sqff = Sqff_cache_data(std::list<std::pair<Polynomial_2,int> >());
        m_cak_2->m_algebraic_kernel_2.square_free_factorize_2_object()
          (p,std::back_inserter(*sqff));
      }
      std::copy(sqff->begin(),sqff->end(),oit);
    }
  };
  
  Square_free_factorize_2 square_free_factorize_2_object() const { 
    return Square_free_factorize_2(this); 
  }


};


} //namespace CGAL



#endif // CGAL_CHACHED_ALGEBRAIC_KERNEL_2_H
