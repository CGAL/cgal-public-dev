// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_RATIONAL_POINT_ABOVE_OR_BELOW_ARC_H
#define CGAL_RATIONAL_POINT_ABOVE_OR_BELOW_ARC_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/VDOL_3/get_boundary_in_range.h>

namespace CGAL {

namespace VDOL_3
{

//! Construct a rational point above or below an algebraic curve x-monotone arc.
/*! Construct a rational point above or below an algebraic curve x-monotone arc.
  "Above" is defined to be the area to the area to the left of the curve when
  traveling on it from a lexicographic smaller point to a lexicographic larger
  point. The function handles vertical arcs.
  \param arc X-monotone arc
  \param above True if the function should construct an above point. False
  if the function should construct a below point.
  \return A rational pair (representing x and y respectively) which forms a 
  rational point above/below the given x-monotone arc.
  \todo check why can't we use construct_interior_vertex
*/

// Within the overlay of the given guards and the supporting curve of the arc
// this function computes a rational point in the cell above the given arc. 
// Precondition: 
//   the arc may intesect the guards only at its endpoints. 
//   or is supported by a guard   

// template <class SVCET_3, class InputIterator>
// std::pair< typename SVCET_3::FT, typename SVCET_3::FT >
// rational_point_above_or_below_arc(
//     const SVCET_3* svcet_3, 
//     const typename SVCET_3::X_monotone_curve_2& arc, 
//     InputIterator guards_begin, InputIterator guards_end, 
//     bool above)
// {
//   typedef typename SVCET_3::FT FT; 
//   std::pair<FT,FT> key = _rational_point_above_or_below_arc(svcet_3,arc,guards_begin, guards_end,above); 
//   typename SVCET_3::Cache_container::Rational_point_above_or_below_arc_cache_data& data = 
//     svcet_3->caches_ptr()->m_rational_point_above_or_below_arc_cache[key];
//   std::cout << data.is_initialized() << std::flush; 
//   data = key; 
//   return key;
// }
         
template <class SVCET_3, class InputIterator>
std::pair< typename SVCET_3::FT, typename SVCET_3::FT >
rational_point_above_or_below_arc(
    const SVCET_3* svcet_3, 
    const typename SVCET_3::X_monotone_curve_2& arc, 
    InputIterator guards_begin, InputIterator guards_end, 
    bool above)
{
  CGAL_SNAP_SVCET_3_TYPEDEFS;

  typedef typename SVCET_3::FT FT; 
  typedef std::pair<FT,FT > result_type;
  
  typedef typename SVCET_3::Algebraic_kernel_d_2 Algebraic_kernel_d_2;
  typedef typename SVCET_3::Algebraic_kernel_d_1 Algebraic_kernel_d_1;
  typedef typename Algebraic_kernel_d_2::Algebraic_real_2 Algebraic_real_2;
  typedef typename Algebraic_kernel_d_2::Algebraic_real_1 Algebraic_real_1;
  typedef typename Algebraic_kernel_d_1::Solve_1 Solve_1;


  Algebraic_kernel_d_2 ak = svcet_3->kernel();
//  Solve_1 solve_1 = ak.solve_1_object();

  typename CGAL::Fraction_traits<FT>::Decompose decompose; 

//   std::cout <<"  arc     : " << arc.curve().polynomial_2() << std::endl;
//   std::cout <<"            " << arc << std::endl;
//   std::cout <<"  guards : " << std::endl;
//   std::copy(guards_begin, guards_end, std::ostream_iterator<Poly_int_2>(std::cout,"\n"));
//   std::cout << std::flush;
    

  if (!arc.is_vertical()){
#if 0
    // currently this code does not work since AK_2::Isolate_2(a,f,g) expects 
    // f(a) != 0 and g(a) != 0
    assert(0); 
    result_type result;
    
    result.first = get_boundary_in_x_range_interior(svcet_3,arc);
    Integer x_num, x_den; 
    decompose(result.first,x_num,x_den);    
    Poly_int_2 slice(Poly_int_1(x_den,x_num));

    typename SVCET_3::Point_2 p(Coordinate_1(result.first),arc.curve(),arc.arcno());    
    typename Algebraic_kernel_d_2::Isolate_2::result_type box; 
    
    box = svcet_3->kernel().isolate_2_object()(p.xy(),arc.curve().polynomial_2(),slice);
    result.second = above ? box[3] : box[2];
    
    for(InputIterator it = guards_begin; it != guards_end; it++){
      box = svcet_3->kernel().isolate_2_object()(p.xy(),*it,slice); 
      result.second = above ? 
        (std::max)(box[3],result.second) : 
        (std::min)(box[2],result.second);
    }
    return result; 

#else
    FT x_rat = get_boundary_in_x_range_interior(svcet_3,arc);
    
    Integer x_num, x_den; 
    decompose(x_rat,x_num,x_den);    
    
    typename SVCET_3::Poly_int_1 slice = 
      CGAL::canonicalize(
        CGAL::evaluate_homogeneous(
             svcet_3->swap_2(arc.curve().polynomial_2()),x_num,x_den));
    
    std::vector<Algebraic_real_1> roots;
    // solve_1(slice,std::back_inserter(roots),false); 
    svcet_3->solve_1(slice,std::back_inserter(roots));
    std::sort(roots.begin(),roots.end());
    Algebraic_real_1 y_on_arc = roots[arc.arcno()];
 
    // put all roots in one container
    for(InputIterator it = guards_begin; it != guards_end; it++){
      typename SVCET_3::Poly_int_1 guard_slice = 
        CGAL::canonicalize(
            CGAL::evaluate_homogeneous(
                svcet_3->swap_2(*it),x_num,x_den));
      // solve_1(guard_slice,std::back_inserter(roots),false);    
      svcet_3->solve_1(guard_slice,std::back_inserter(roots));
    } 
  
    Algebraic_real_1 other;
    if(above){
      // we have to find a smallest number larger than y_on_arc
      other = *(std::max_element(roots.begin(),roots.end()));
      
      for(unsigned int i = 0; i < roots.size(); i++){
        if(roots[i] > y_on_arc ){
          other = min(roots[i],other);
          CGAL_postcondition(other > y_on_arc ); 
        }
      }
      if(other==y_on_arc){
        return std::make_pair(x_rat,FT(ak.upper_boundary_1_object()(y_on_arc)) + FT(1));
      }else{
        CGAL_postcondition(other > y_on_arc ); 
        return std::make_pair(x_rat,FT(ak.boundary_between_1_object()(y_on_arc,other)));
      }
    }else{
      // we have to find a largest number smaller than y_on_arc
      other = *(std::min_element(roots.begin(),roots.end()));
      
      for(unsigned int i = 0; i < roots.size(); i++){
        if(roots[i] < y_on_arc )
          other = (max)(roots[i],other);
      }
      if(other==y_on_arc){
        // arc is the lower most arc
        return std::make_pair(x_rat,FT(ak.lower_boundary_1_object()(y_on_arc)) - FT(1));
      }else{
        CGAL_postcondition(other < y_on_arc ); 
        return std::make_pair(x_rat,FT(ak.boundary_between_1_object()(y_on_arc,other)));
      }
    }
#endif
  }else{
    FT y_rat = get_boundary_in_y_range_interior(svcet_3,arc);
  
    
    Integer y_num, y_den; 
    decompose(y_rat,y_num,y_den);    
    
    
    typename SVCET_3::Poly_int_1 slice = 
      CGAL::canonicalize(
        CGAL::evaluate_homogeneous(arc.curve().polynomial_2(),y_num,y_den));
    
    

    std::vector<Algebraic_real_1> roots;
    //solve_1(slice,std::back_inserter(roots),false); 
    svcet_3->solve_1(slice,std::back_inserter(roots));

    CGAL_assertion(std::find(roots.begin(),roots.end(),arc.x())!=roots.end());

    // put all roots in one container
    for(InputIterator it = guards_begin; it != guards_end; it++){
      typename SVCET_3::Poly_int_1 guard_slice = 
        CGAL::canonicalize(
            CGAL::evaluate_homogeneous(*it,y_num,y_den));
      // solve_1(guard_slice,std::back_inserter(roots),false);    
      svcet_3->solve_1(guard_slice,std::back_inserter(roots));
    } 
  
    Algebraic_real_1 other;
    // this is all in vertical direction 
    // hence above means to the left, which is smaller
    if(above){
      other = *(std::min_element(roots.begin(),roots.end()));
      for(unsigned int i = 0; i < roots.size(); i++){
        if(roots[i] < arc.x() )
          other = (max)(roots[i],other);
      }
      CGAL_assertion(other <= arc.x());
      if(other == arc.x()){
        return std::make_pair(FT(ak.lower_boundary_1_object()(arc.x())) - FT(1),y_rat);
      }else{
        CGAL_postcondition(other < arc.x() ); 
        return std::make_pair(FT(ak.boundary_between_1_object()(other,arc.x())),y_rat);
      }
    }else{
      other = *(std::max_element(roots.begin(),roots.end()));
      for(unsigned int i = 0; i < roots.size(); i++){
        if(roots[i] > arc.x() )
          other = (min)(roots[i],other);
      }
      CGAL_assertion(other >= arc.x());
      if(other == arc.x()){
        return std::make_pair(FT(ak.upper_boundary_1_object()(arc.x())) + FT(1),y_rat);
      }else{
        CGAL_postcondition(other > arc.x() ); 
        return std::make_pair(FT(ak.boundary_between_1_object()(arc.x(),other)),y_rat);
      }
    }
  }  
}

} // namespace VDOL_3
} //namespace CGAL

#endif // CGAL_RATIONAL_POINT_ABOVE_OR_BELOW_ARC_H
