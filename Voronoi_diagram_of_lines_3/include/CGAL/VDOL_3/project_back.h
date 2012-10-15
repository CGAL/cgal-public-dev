
#ifndef CGAL_VDOL_PROJECT_BACK_H
#define CGAL_VDOL_PROJECT_BACK_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/envelope_3.h>
#include <CGAL/VDOL_3/vol_tools.h>
#include <CGAL/VDOL_3/svcet_3_state_dependent_functions.h>
#include <CGAL/VDOL_3/CGAL_SNAP_SVCET_3_TYPEDEFS.h>
#include <CGAL/intersections.h>
#include <CGAL/squared_distance_3.h>


namespace CGAL {
namespace VDOL_3 { 
      
template <class SVCET_3>
typename SVCET_3::Point_3 
project_back(
    const SVCET_3& svcet_3, 
    const typename SVCET_3::Point_2& point_2, 
    const typename SVCET_3::Poly_int_3& tbs) {
  
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  
  std::pair<double,double> approx = point_2.xy().to_double();
  FT u(approx.first);
  FT v(approx.second);
  FT r; 

  typename CGAL::Coercion_traits<Poly_int_3,Poly_rat_3>::Cast cast; 
  Poly_rat_3 poly_3 = cast(CGAL::swap(tbs,0,2));
  Poly_rat_1 poly_1 = CGAL::canonicalize(CGAL::evaluate(CGAL::evaluate(poly_3,Poly_rat_2(u)),Poly_rat_1(v)));
  CGAL_postcondition(CGAL::is_zero(poly_1) || poly_1.lcoeff()==1);

  switch(CGAL::degree(poly_1)){
  case 0: 
    // return Point_3(1000,1000,1000);
    r = 0; 
    break; 
  case 1: 
    r = -poly_1[0]/poly_1[1];
    // bisector is just a plane and the point is on the wrong side
    // or the point is on the infinite line at infinity. 
    if(r < 0) r = 0; // return Point_3(1000,1000,1000);
    break; 
  default: // = 2 
    CGAL_postcondition(CGAL::degree(poly_1) == 2); 
    FT p = poly_1[1];
    FT q = poly_1[0];
    FT discr = p*p/4 - q;
    if(discr < 0) r = 10; //  Point_3(1000,1000,1000);
                                                                        \
    double sdiscr = CGAL::sqrt(CGAL::to_double(CGAL::abs(discr)));
    r = FT(-p/2 + sdiscr);
     
  }

  typedef typename SVCET_3::Linear_kernel Linear_kernel; 
  Linear_kernel const& lk = svcet_3.linear_kernel(); 

  Vector_3 v1 = lk.construct_scaled_vector_3_object()(svcet_3.v1(),u);
  Vector_3 v2 = lk.construct_scaled_vector_3_object()(svcet_3.v2(),v*r);
  Vector_3 v3 = lk.construct_scaled_vector_3_object()(svcet_3.v3(),r);
  Point_3  ep = svcet_3.origin();

  ep = lk.construct_translated_point_3_object()(ep,v1);
  ep = lk.construct_translated_point_3_object()(ep,v2);
  ep = lk.construct_translated_point_3_object()(ep,v3);

  return ep; 
}


template <class SVCET_3, class OutputIterator>
OutputIterator 
project_back(
    const SVCET_3& svcet_3, 
    const typename Envelope_diagram_2<SVCET_3>::Vertex& v, 
    const Approximation_info& info,
    OutputIterator oit)
{
  if(!v.point().is_finite()) return oit; 
  
  CGAL_SNAP_SVCET_3_TYPEDEFS;

  typedef CGAL::Envelope_diagram_2<SVCET_3> Envelope_diagram_2; 
  typedef typename Envelope_diagram_2::Surface_const_iterator SIT;
  typedef typename Envelope_diagram_2::Xy_monotone_surface_3 Surface; 
    
  CGAL_precondition(v.number_of_surfaces()>0);    
  //std::cout << " \n PROCESS VERTEX \n "  << std::endl; 

  SIT singular_surface = v.surfaces_end();
  for(SIT sit = v.surfaces_begin(); sit != v.surfaces_end(); sit++){
    if(is_singular_line_at(&svcet_3,v.point(),*sit))
      singular_surface = sit;  
  }
 
  if(singular_surface == v.surfaces_end()){
    // point at infitiy or normal 
    if(has_infinit_distance_at(&svcet_3,v.point(),v.surface()))
      return oit; 
    else 
      return *oit++ = project_back(svcet_3,v.point(),v.surface().transformed_bisector());
  }

  
  SIT min_surface = v.surfaces_end(); // that is not a singular line 
  for(SIT sit = v.surfaces_begin(); sit != v.surfaces_end(); sit++){ 
    if(!is_singular_line_at(&svcet_3,v.point(),*sit)){
      if( min_surface == v.surfaces_end()){
        min_surface = sit;
      }else{
        if( svcet_3.compare_z_at_xy_3_object()(v.point(),*sit,*min_surface) == CGAL::SMALLER) 
          min_surface = sit;
      }
    }
  }
  
  // three options 
  // min_surface != v.surfaces_end() && min_surface has zero distance -> Vertex on base line -> p_max = p_min 
  // min_surface != v.surfaces_end() && min_surface is normal -> EDGE -> pmax on min surface != p_min on base line 
  // min_surface == v.surfaces_end() || has_infinit_distance_at -> RAY  -> p_max defined by clipping sphere 

  bool is_ray  = false; 
  if( min_surface != v.surfaces_end()){
    if(has_zero_distance_at(&svcet_3,v.point(),*min_surface)){
      // is a vertex on the base line 
      return *oit++ = project_back(svcet_3,v.point(),min_surface->transformed_bisector());
    }else{
      if(has_infinit_distance_at(&svcet_3,v.point(),*min_surface)){
        is_ray = true; 
      } // else is edge 
    }
  }else{
    is_ray = true; 
  }
 
  Point_3 p_min,p_max; 
  if (is_ray){
    //std::cout << "IS RAY" << std::endl; 
    p_min =  project_back(svcet_3,v.point(),singular_surface->transformed_bisector()); 
    Poly_int_3 x = shift(Poly_int_3(1),1,0);
    Poly_int_3 y = shift(Poly_int_3(1),1,1);
    Poly_int_3 z = shift(Poly_int_3(1),1,2);
    Poly_int_3 csphere = x*x+y*y+z*z - info.sradius_clipping_sphere; 
    Poly_int_3 tcsphere = VDOL_3::construct_transformed_bisector_3(&svcet_3,csphere);   
    p_max =   project_back(svcet_3,v.point(),tcsphere);    
  }else{ // is_edge
    //std::cout << "IS EDGE" << std::endl; 
    CGAL_precondition( min_surface != v.surfaces_end());
    p_min =  project_back(svcet_3,v.point(),singular_surface->transformed_bisector()); 
    p_max =  project_back(svcet_3,v.point(),min_surface->transformed_bisector()); 
  }

  //std::cout << "p_min: " << CGAL::to_double(p_min.x()) <<" " << CGAL::to_double(p_min.y()) <<" " << CGAL::to_double(p_min.z()) << "\t " << p_min  <<" " << std::endl;
  //std::cout << "p_max: " << CGAL::to_double(p_max.x()) <<" " << CGAL::to_double(p_max.y()) <<" " << CGAL::to_double(p_max.z()) << "\t " << p_max  <<" " << std::endl;

  double distance = CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(p_max,p_min)));
  double delta    = CGAL::sqrt(info.sdistance_generation);
  double rdelta   = delta/distance; 

  //std::cout << distance << "\t" << delta << "\t"<< rdelta << std::endl; 
  
  *oit++ = p_min; 
  *oit++ = p_max; 
  if( distance * distance <= info.sradius_clipping_sphere + 0.1 ){
    for(double t = rdelta; t < 1; t+= rdelta ){
      Point_3 p(p_max.x()*(1-t)+p_min.x()*t,p_max.y()*(1-t)+p_min.y()*t,p_max.z()*(1-t) + p_min.z()*t);
      *oit++=(p);
    } 
  }
  return oit; 
}


template <class SVCET_3, class OutputIterator>
OutputIterator split_vertical_arc(
    const SVCET_3& svcet_3, 
    const typename SVCET_3::Coordinate_1 x,
    const typename SVCET_3::Rational& y_min,
    const typename SVCET_3::Rational& y_max, 
    const typename SVCET_3::Poly_int_3& tbs,
    const Approximation_info& info, 
    OutputIterator oit)   
{
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  Rational y_mid((y_min+y_max)/2);

  Point_3 p_min_3 = project_back(svcet_3,construct_point_2(&svcet_3,x,y_min),tbs);
  Point_3 p_mid_3 = project_back(svcet_3,construct_point_2(&svcet_3,x,y_mid),tbs);
  Point_3 p_max_3 = project_back(svcet_3,construct_point_2(&svcet_3,x,y_max),tbs);
    
  if( CGAL::squared_distance(p_min_3,p_max_3) < info.sdistance_generation)
    return oit++ =  p_min_3;
  if( to_double(y_max)-to_double(y_min) < 0.00001 ) 
    return oit++ =  p_min_3;

  bool is_visible_min = (slength(p_min_3) <= info.sradius_bounding_sphere );
  bool is_visible_mid = (slength(p_mid_3) <= info.sradius_bounding_sphere );
  bool is_visible_max = (slength(p_max_3) <= info.sradius_bounding_sphere );
    
  if(is_visible_min || is_visible_mid) 
    split_vertical_arc(svcet_3,x,y_min,y_mid,tbs,info,oit);
  if(is_visible_mid || is_visible_max)
    split_vertical_arc(svcet_3,x,y_mid,y_max,tbs,info,oit);
  return oit;
}




template <class SVCET_3, class OutputIterator>
OutputIterator split_arc(
    const SVCET_3& svcet_3, 
    const typename SVCET_3::X_monotone_curve_2& xcurve, 
    const typename SVCET_3::Poly_int_3& tbs,
    const typename SVCET_3::Coordinate_1& x_min,
    const typename SVCET_3::Coordinate_1& x_max, 
    const typename SVCET_3::Point_3& p_min_3, 
    const typename SVCET_3::Point_3& p_max_3,
    const Approximation_info& info, 
    OutputIterator oit) {
  CGAL_SNAP_SVCET_3_TYPEDEFS;
 
  Coordinate_1 x_mid(svcet_3.curved_kernel_2().kernel().bound_between_1_object()(x_min,x_max));
  Point_2 p_mid_2 (x_mid,xcurve.curve(),xcurve.arcno());
  Point_3 p_mid_3 = project_back(svcet_3,p_mid_2,tbs);
  
  if(CGAL::squared_distance(p_min_3,p_max_3) < info.sdistance_generation  || to_double(x_max)-to_double(x_min) < 0.0000001) {
    return *oit++ = p_min_3;
  }
    
  bool is_visible_min = (slength(p_min_3) <= info.sradius_bounding_sphere );
  bool is_visible_mid = (slength(p_mid_3) <= info.sradius_bounding_sphere );
  bool is_visible_max = (slength(p_max_3) <= info.sradius_bounding_sphere );
    
  if(is_visible_min || is_visible_mid)
    oit = split_arc(svcet_3,xcurve,tbs,x_min,x_mid,p_min_3,p_mid_3,info,oit);
  if(is_visible_mid || is_visible_max)
    oit = split_arc(svcet_3,xcurve,tbs,x_mid,x_max,p_mid_3,p_max_3,info,oit);
  return oit;
}
 


template <class SVCET_3, class OutputIterator>
OutputIterator
project_back(
    const SVCET_3& svcet_3,  
    typename SVCET_3::X_monotone_curve_2 xcurve,
    typename SVCET_3::Poly_int_3 tbs,
    const Approximation_info& info, 
    OutputIterator oit)
{

  //std::cout << "." << std::endl;   
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  
// typename SVCET_3::X_monotone_curve_2 xcurve = e.curve();
//   Poly_int_3 tbs = e.surface().transformed_bisector();
  
  typename Linear_kernel::Compute_squared_length_3 slength 
    = svcet_3.linear_kernel().compute_squared_length_3_object() ;
  typename Linear_kernel::Construct_vector_3 construct_v 
    = svcet_3.linear_kernel().construct_vector_3_object();

  if(xcurve.is_vertical()){
    //std::cout <<" split vertical arc "<< std::endl; 
    Coordinate_1 x = xcurve.curve_end_x(ARR_MAX_END);
    Rational y_min,y_max;
    if(xcurve.location(ARR_MIN_END) == ARR_BOTTOM_BOUNDARY){
      y_min = -100000; 
    }else{
      y_min = xcurve.curve_end(ARR_MIN_END).xy().y().upper();
    }
  
    if(xcurve.location(ARR_MAX_END) == ARR_TOP_BOUNDARY){
      y_max =  100000; 
    }else{
      y_max = xcurve.curve_end(ARR_MAX_END).xy().y().lower();
    }     
      
    oit++ = project_back(svcet_3,construct_point_2(&svcet_3,x,y_min),tbs);      
    oit++ = project_back(svcet_3,construct_point_2(&svcet_3,x,y_max),tbs);

    if( CGAL::sign(y_min) * CGAL::sign(y_max) == CGAL::NEGATIVE ){
      oit = split_vertical_arc(svcet_3,x,y_min,Rational(0), tbs,info,oit);
      oit = split_vertical_arc(svcet_3,x,Rational(0),y_max,tbs,info, oit);
    }else{      
      oit = split_vertical_arc(svcet_3,x,y_min,y_max,tbs,info, oit);
    }
    //std::cout <<" split vertical arc done "<< std::endl; 
    return oit;
  }

    
  Coordinate_1 x_min,x_max;
  Point_3 p_min_3, p_max_3; 
    
  Arr_parameter_space min_loc = xcurve.location(ARR_MIN_END);
  if(min_loc == ARR_INTERIOR ){
    x_min =  xcurve.curve_end_x(ARR_MIN_END);
    //p_min_3 = project_back(svcet_3,xcurve.curve_end(ARR_MIN_END),tbs);
  }else{
    if( min_loc == ARR_BOTTOM_BOUNDARY || min_loc == ARR_TOP_BOUNDARY){
      x_min = xcurve.curve_end_x(ARR_MIN_END);
      x_min = Coordinate_1(svcet_3.curved_kernel_2().kernel().approximate_relative_1_object()(x_min,53).second+0.000001);
    }else{
      x_min = Coordinate_1(-110);
    }
    //p_min_3 = project_back(svcet_3,Point_2(x_min,xcurve.curve(),xcurve.arcno()),tbs);
  }

  Arr_parameter_space max_loc = xcurve.location(ARR_MAX_END);
  if(max_loc == ARR_INTERIOR ){
    x_max =  xcurve.curve_end_x(ARR_MAX_END);
    //p_max_3 = project_back(svcet_3,xcurve.curve_end(ARR_MAX_END),tbs);
  }else{
    if( max_loc == ARR_BOTTOM_BOUNDARY || max_loc == ARR_TOP_BOUNDARY){
      x_max = xcurve.curve_end_x(ARR_MAX_END);
      x_max = Coordinate_1(svcet_3.curved_kernel_2().kernel().approximate_relative_1_object()(x_max,53).first-0.000001);
    }else{  
      x_max = Coordinate_1(110);
    }
    //p_max_3 = project_back(svcet_3,Point_2(x_max,xcurve.curve(),xcurve.arcno()),tbs); 
  }
  if(x_min > x_max) {
//      std::cerr << "I had to skip an edge because it seems too small" << std::endl;      
//      std::cerr << to_double(x_min) << " " << to_double(x_max) << std::endl;
    return oit; // in case an edge is too small 
  }
    
    
  if(!xcurve.is_in_x_range_interior(x_min)){
    p_min_3 = project_back(svcet_3,xcurve.curve_end(ARR_MIN_END),tbs);
  }else{
    p_min_3 = project_back(svcet_3,Point_2(x_min,xcurve.curve(),xcurve.arcno()),tbs); 
  }
    
  if(!xcurve.is_in_x_range_interior(x_max)){
    p_max_3 = project_back(svcet_3,xcurve.curve_end(ARR_MAX_END),tbs);
  }else{
    p_max_3 = project_back(svcet_3,Point_2(x_max,xcurve.curve(),xcurve.arcno()),tbs); 
  }
 
  *oit++ = p_min_3;
  *oit++ = p_max_3; 

  
  if( CGAL::sign(x_min) * CGAL::sign(x_max) == CGAL::NEGATIVE ){      
    Coordinate_1 x_mid(0);
    Point_3 p_mid_3(project_back(svcet_3,Point_2(x_mid,xcurve.curve(),xcurve.arcno()),tbs)); 
    oit =  split_arc(svcet_3,xcurve,tbs,x_min,x_mid,p_min_3,p_mid_3,info,oit);
    oit =  split_arc(svcet_3,xcurve,tbs,x_mid,x_max,p_mid_3,p_max_3,info,oit);
  }else{
    oit =  split_arc(svcet_3,xcurve,tbs,x_min,x_max,p_min_3,p_max_3,info,oit);
  }
  return oit; 
}

template <class SVCET_3, class OutputIterator>
OutputIterator
project_back(
    const SVCET_3& svcet_3, 
    const typename Envelope_diagram_2<SVCET_3>::Halfedge& e, 
    const Approximation_info& info, 
    OutputIterator oit)
{
  return project_back(svcet_3,e.curve(),e.surface().transformed_bisector(),info,oit); 
} 


// template <class SVCET_3, class OutputIterator>
// OutputIterator
// project_back(
//     const SVCET_3& svcet_3, 
//     const typename Envelope_diagram_2<SVCET_3>::Face& f, 
//     const Approximation_info& info,
//     OutputIterator oit){
//   return oit; 
// }





} // namespace VDOL_3
} // namespace CGAL 

#endif // CGAL_VDOL_PROJECT_BACK_H
