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

#ifndef CGAL_VOL_TOOLS_H
#define CGAL_VOL_TOOLS_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Arr_enums.h>
#include <boost/optional.hpp>
#include <vector>

namespace CGAL {

namespace VDOL_3
{


template < typename Halfedge>
bool is_vertical(Halfedge* he){
  typedef typename Halfedge::Vertex Vertex; 

  CGAL_assertion(he->has_null_curve());

  const Vertex    *v1 = he->vertex();
  const Vertex    *v2 = he->opposite()->vertex();
  
  if( v1->parameter_space_in_x() == ARR_LEFT_BOUNDARY && 
      v2->parameter_space_in_x() == ARR_LEFT_BOUNDARY ) 
    return true;
  if( v1->parameter_space_in_x() == ARR_RIGHT_BOUNDARY && 
      v2->parameter_space_in_x() == ARR_RIGHT_BOUNDARY ) 
    return true;
  return false;
}
template <typename Vertex>
int degree(Vertex  *v){
  if(v->halfedge() == NULL) return 0; 
  
  typename Vertex::Halfedge  *he = v->halfedge();
  he = he->next()->opposite();
  int result = 1; 

  while(he != v->halfedge()){
    he = he->next()->opposite();
    result++;
  }
  return result;
}



template <typename SVCET_3, typename OutputIterator>
OutputIterator construct_x_monotone_curves_2(
    const SVCET_3* svcet_3,
    const typename SVCET_3::Poly_int_2 poly,
    OutputIterator oit){

  // TODO: consider caching, the traits is present.  
  
  typename SVCET_3::Curve_2 curve =  
    svcet_3->kernel().construct_curve_2_object()(poly);
  
  typename SVCET_3::Obj_list x_mono_objects;
  svcet_3->make_x_monotone_2_object()
    (curve, std::back_inserter(x_mono_objects));
 
  typedef typename SVCET_3::X_monotone_curve_2 X_monotone_curve_2;
  std::list<X_monotone_curve_2> xcurves; 
  for (typename SVCET_3::Obj_list::iterator it = x_mono_objects.begin();
       it != x_mono_objects.end(); ++it){
    typename SVCET_3::X_monotone_curve_2 xcurve;
    CGAL_assertion(assign(xcurve, *it));
    assign(xcurve, *it);
    *oit++ = xcurve;
  }
  return oit; 
}



/*
  The function construct_bisector_3 constructs the bisector of two lines in R3
  In general, the resulting polynomial represents a quadric surface. 
  The code was generated with the help of the subsequent maple code:


  with(LinearAlgebra); with(codegen, C);
  p1 := Vector([AA, BB, CC]); v1 := Vector([DD, EE, FF]); 
  p2 := Vector([GG, HH, II]); v2 := Vector([JJ, KK, LL]);
  X := Vector([x,y,z]):
  D1 := DotProduct( v1 &x (X - p1) , v1 &x (X - p1),conjugate=false) :
  D2 := DotProduct( v2 &x (X - p2) , v2 &x (X - p2),conjugate=false) :
  Q2 := expand(D1*DotProduct(v2, v2, conjugate = false)
        -D2*DotProduct(v1, v1, conjugate = false)):
  a := coeff(coeff(coeff(Q2,x,2),y,0),z,0):
  b := coeff(coeff(coeff(Q2,x,0),y,2),z,0):
  c := coeff(coeff(coeff(Q2,x,0),y,0),z,2):

  d := coeff(coeff(coeff(Q2,x,1),y,1),z,0):
  e := coeff(coeff(coeff(Q2,x,1),y,0),z,1):
  f := coeff(coeff(coeff(Q2,x,0),y,1),z,1):

  g := coeff(coeff(coeff(Q2,x,1),y,0),z,0):
  h := coeff(coeff(coeff(Q2,x,0),y,1),z,0):
  i := coeff(coeff(coeff(Q2,x,0),y,0),z,1):

  j := coeff(coeff(coeff(Q2,x,0),y,0),z,0):

  A := array([a,b,c,d,e,f,g,h,i,j]):
  C(A,optimized);
*/
template <typename Kernel>
Polynomial<Polynomial<Polynomial<typename Kernel::FT> > >
construct_bisector_3(
    const typename Kernel::Line_3& l1, 
    const typename Kernel::Line_3& l2){

  // CGAL_precondition(!typename Kernel::Do_intersect_3()(l1,l2));

  typedef typename Kernel::FT FT;

  Kernel ker;

  typename Kernel::Construct_point_on_3 get_point = 
    ker.construct_point_on_3_object();
  typename Kernel::Compute_x_3 compute_x_3 = ker.compute_x_3_object();
  typename Kernel::Compute_y_3 compute_y_3 = ker.compute_y_3_object();
  typename Kernel::Compute_z_3 compute_z_3 = ker.compute_z_3_object();
  typename Kernel::Compute_dx_3 compute_dx_3 = ker.compute_dx_3_object();
  typename Kernel::Compute_dy_3 compute_dy_3 = ker.compute_dy_3_object();
  typename Kernel::Compute_dz_3 compute_dz_3 = ker.compute_dz_3_object();
  typename Kernel::Construct_direction_3 get_direction = 
    ker.construct_direction_3_object();
    
  typename Kernel::Point_3 p1     = get_point(l1,0);
  typename Kernel::Point_3 p2     = get_point(l2,0);
  typename Kernel::Direction_3 v1 = get_direction(l1);
  typename Kernel::Direction_3 v2 = get_direction(l2);

  FT AA = compute_x_3(p1); 
  FT BB = compute_y_3(p1); 
  FT CC = compute_z_3(p1);
  FT DD = compute_dx_3(v1); 
  FT EE = compute_dy_3(v1); 
  FT FF = compute_dz_3(v1);  
  
  FT GG = compute_x_3(p2); 
  FT HH = compute_y_3(p2); 
  FT II = compute_z_3(p2);
  FT JJ = compute_dx_3(v2); 
  FT KK = compute_dy_3(v2); 
  FT LL = compute_dz_3(v2);
    
  FT t1 = FF*FF;
  FT t2 = JJ*JJ;
  FT t3 = t1*t2;
  FT t4 = EE*EE;
  FT t5 = t4*t2;
  FT t6 = KK*KK;
  FT t7 = DD*DD;
  FT t8 = t6*t7;
  FT t9 = LL*LL;
  FT t10 = t9*t7;
  FT t12 = t1*t6;
  FT t13 = t9*t4;
  FT t16 = DD*EE;
  FT t20 = JJ*KK;
  FT t25 = FF*DD;
  FT t29 = LL*JJ;
  FT t34 = EE*FF;
  FT t37 = KK*LL;
  FT t45 = DD*BB;
  FT t48 = t1*AA;
  FT t52 = t4*AA;
  FT t56 = t6*GG;
  FT t59 = t9*GG;
  FT t62 = t25*CC*t9+t45*EE*t2-t48*t9-t48*t6-t48*t2-t52*t9-t52*t2-t52*t6+t56*
    t4+t56*t1+t59*t1+t59*t4;
  FT t78 = JJ*HH;
  FT t86 = t59*t7+t25*CC*t2+t25*CC*t6-t29*II*t1+t45*EE*t6+t45*EE*t9-t29*II*t7-
    t29*II*t4-t78*KK*t7-t78*KK*t4+t56*t7-t78*KK*t1;
  FT t88 = KK*II;
  FT t93 = t2*HH;
  FT t96 = t9*HH;
  FT t100 = t7*BB;
  FT t108 = -t88*LL*t7-t20*GG*t1+t93*t4+t93*t1+t96*t7+t96*t1+t96*t4-t100*t9-
    t100*t2-t100*t6-t20*GG*t7-t88*LL*t4;
  FT t113 = t1*BB;
  FT t120 = EE*CC;
  FT t131 = -t88*LL*t1-t20*GG*t4-t113*t9-t113*t6-t113*t2+t93*t7+t16*AA*t2+t120
    *FF*t6+t16*AA*t6+t16*AA*t9+t120*FF*t2+t120*FF*t9;
  FT t133 = t4*CC;
  FT t137 = t7*CC;
  FT t141 = FF*AA;
  FT t154 = -t133*t6-t133*t2-t133*t9-t137*t2-t137*t6-t137*t9+t141*DD*t2+t141*
    DD*t6+t141*DD*t9+t34*BB*t2+t34*BB*t6+t34*BB*t9;
  FT t155 = t6*II;
  FT t159 = LL*GG;
  FT t172 = t2*II;
  FT t176 = t155*t4-t37*HH*t7-t159*JJ*t7-t159*JJ*t4-t159*JJ*t1+t155*t7+t155*t1
    -t37*HH*t4-t37*HH*t1+t172*t7+t172*t4+t172*t1;
  FT t178 = GG*GG;
  FT t179 = t9*t178;
  FT t181 = FF*BB;
  FT t185 = CC*CC;
  FT t186 = t4*t185;
  FT t189 = BB*BB;
  FT t190 = t1*t189;
  FT t194 = AA*AA;
  FT t195 = t1*t194;
  FT t197 = HH*HH;
  FT t198 = t2*t197;
  FT t201 = t7*t185;
  FT t204 = t7*t189;
  FT t206 = -t179*t7-FT(2)*t120*t181*t9+t186*t2+t186*t6+t190*t2+t190*t6+t190*t9+
    t195*t2-t198*t7+t195*t9+t201*t2+t201*t9+t204*t2;
  FT t209 = t4*t194;
  FT t213 = II*II;
  FT t214 = t6*t213;
  FT t217 = t9*t197;
  FT t220 = DD*CC;
  FT t228 = t204*t6+t204*t9+t209*t2+t209*t6+t209*t9-t214*t4-t214*t1-t217*t4-
    t217*t1-FT(2)*t141*t220*t2-t179*t4-t217*t7-t214*t7-t179*t1;
  FT t234 = LL*HH;
  FT t241 = JJ*II;
  FT t252 = KK*GG;
  FT t259 = -t141*t220*t6-t141*t220*t9+t88*t234*t7+t88*t234*t4+t88*t234*t1+
    t159*t241*t7-t120*t181*t2-t120*t181*t6+t159*t241*t4+t159*t241*t1+t78*t252*t7+
    t78*t252*t4+t78*t252*t1;
  FT t260 = EE*AA;
  FT t270 = t2*t213;
  FT t276 = t6*t178;
  FT t283 = -FT(2)*t45*t260*t2-FT(2)*t45*t260*t6-FT(2)*t45*t260*t9-t270*t7-t270*t4-
    t270*t1-t198*t4-t198*t1-t276*t7-t276*t1+t201*t6+t195*t6-t276*t4+t186*t9;

  FT a = t3+t5-t8-t10;
  FT b = t12+t8-t13-t5;
  FT c = t13+t10-t12-t3;
  FT d = -FT(2)*t16*t6-FT(2)*t16*t9-FT(2)*t16*t2+FT(2)*t20*t1+FT(2)*t20*t7+FT(2)*t20*t4;
  FT e  = -FT(2)*t25*t6-FT(2)*t25*t2-FT(2)*t25*t9+FT(2)*t29*t7+FT(2)*t29*t4+FT(2)*t29*t1;
  FT f = -FT(2)*t34*t9-FT(2)*t34*t6+FT(2)*t37*t7+FT(2)*t37*t4+FT(2)*t37*t1-FT(2)*t34*t2;
  FT g = FT(2)*t62+FT(2)*t86;
  FT h = FT(2)*t108+FT(2)*t131;
  FT i = FT(2)*t154+FT(2)*t176;
  FT j = t206+t228+FT(2)*t259+t283;

  typedef CGAL::Polynomial<FT> Polynomial_1;
  typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;
  typedef CGAL::Polynomial<Polynomial_2> Polynomial_3;
  typedef CGAL::Polynomial_traits_d<Polynomial_3> PT;

  typedef std::pair<CGAL::Exponent_vector,FT> Monomial;
  std::vector<Monomial> monomials;
  CGAL::Exponent_vector ev(0,0,0);

  ev[0]=2; ev[1]=0; ev[2]=0; monomials.push_back(Monomial(ev,a));
  ev[0]=0; ev[1]=2; ev[2]=0; monomials.push_back(Monomial(ev,b));
  ev[0]=0; ev[1]=0; ev[2]=2; monomials.push_back(Monomial(ev,c)); 

  ev[0]=1; ev[1]=1; ev[2]=0; monomials.push_back(Monomial(ev,d));
  ev[0]=1; ev[1]=0; ev[2]=1; monomials.push_back(Monomial(ev,e));
  ev[0]=0; ev[1]=1; ev[2]=1; monomials.push_back(Monomial(ev,f));

  ev[0]=1; ev[1]=0; ev[2]=0; monomials.push_back(Monomial(ev,g));
  ev[0]=0; ev[1]=1; ev[2]=0; monomials.push_back(Monomial(ev,h));
  ev[0]=0; ev[1]=0; ev[2]=1; monomials.push_back(Monomial(ev,i));

  ev[0]=0; ev[1]=0; ev[2]=0; monomials.push_back(Monomial(ev,j));
  typename PT::Construct_polynomial construct_polynomial;
  Polynomial_3 result = construct_polynomial(monomials.begin(),
      monomials.end());
    
#ifndef CGAL_NDEBUG 

  // postcond: 
  // if the two lines do not intersect, the do not intersect the constructed 
  // bisector. 
  // if they intersect, the bisector consists of two intersection planes and 
  // the lines intersect the planes at their common line (double root)
    
  // typedef typename PT::template Rebind<FT,1>::Other PT_1;
  // typedef typename PT_1::Polynomial_d Polynomial_1;

  CGAL::set_pretty_mode(std::cout);
  std::vector<Polynomial_1> line_param_1;
  line_param_1.push_back(Polynomial_1(AA,DD)); 
  line_param_1.push_back(Polynomial_1(BB,EE));
  line_param_1.push_back(Polynomial_1(CC,FF));
  //for(unsigned int i =0; i < line_param_1.size(); i++){
  //  std::cout << line_param_1[i] <<std::endl;
  //}
  
  std::vector<Polynomial_1> line_param_2;
  line_param_2.push_back(Polynomial_1(GG,JJ)); 
  line_param_2.push_back(Polynomial_1(HH,KK));
  line_param_2.push_back(Polynomial_1(II,LL));
 
  typename CGAL::Polynomial_traits_d<Polynomial_1>::Get_coefficient 
    get_coefficient; 
  Polynomial_1 res1 = 
    typename PT::Substitute()(result,line_param_1.begin(),line_param_1.end());
  // discr = b^2 - 4ac
  FT discr1 = 
    get_coefficient(res1,1)*get_coefficient(res1,1) - 
    FT(4)* get_coefficient(res1,0)*get_coefficient(res1,2); 
  Polynomial_1 res2 = 
    typename PT::Substitute()(result,line_param_2.begin(),line_param_2.end());
  FT discr2 = 
    get_coefficient(res2,1)*get_coefficient(res2,1) - 
    FT(4)* get_coefficient(res2,0)*get_coefficient(res2,2); 
  
  CGAL_postcondition(discr1 <= 0); // won't work for modular arithmetic 
  CGAL_postcondition(discr2 <= 0); // won't work for modular arithmetic 
     
#endif    

  return result;
}


// Given a Point location of a minimization diagram,
// this function computes all surface above the given point via the given output 
// iterator note that the list of surfaces may be empty 
// 
//  Moreover, in case of the vertex things may be more dificult since the 
//  vertex may represent the singular line of a bisector that is induced by 
//  a line that intersects the base line. I this case the minimization diagramm 
//  just stores the intersecting lines which is correct. 
//  However, for the walk of the point location we are actually interested in the 
//  lines that induce a bisector that cuts the singular line. Thus we have to report 
//  lines that induce bisectors in the neighbourhood of the vertex. 
template < class Point_location , class OutputIterator> 
OutputIterator
compute_surfaces(
    const Point_location& pl, 
    const typename Point_location::Point_2& p,
    OutputIterator oit
){
  typedef typename Point_location::Arrangement_2 Envelope_diagram_2; 
  typename Envelope_diagram_2::Face_const_handle       fh;
  typename Envelope_diagram_2::Halfedge_const_handle   eh;
  typename Envelope_diagram_2::Vertex_const_handle     vh;
  typedef typename Envelope_diagram_2::Xy_monotone_surface_3 Surface; 
  typedef typename Envelope_diagram_2::Halfedge_around_vertex_const_circulator 
    Halfedge_around_vertex_const_circulator; 
  
  CGAL::Object obj(pl.locate(p));
  
  if (CGAL::assign (fh, obj)) { 
    if(fh->number_of_surfaces() >0){
      return std::copy(fh->surfaces_begin(),fh->surfaces_end(),oit); 
    }
    return oit; 
  }else if (CGAL::assign (eh, obj)) { 
    if(eh->number_of_surfaces() >0){
      return std::copy(eh->surfaces_begin(),eh->surfaces_end(),oit); 
    }
    return oit; 
  }
  else if (CGAL::assign (vh, obj)){
    std::set<Surface> surfaces; 
    // reporting all surface in the neigborhood.
    std::copy(vh->surfaces_begin(),vh->surfaces_end(),CGAL::inserter(surfaces)); 
    // the minimization diagram should not contain isolated vertices
    CGAL_precondition(!vh->isolated_vertex());
    
    // collect all surfaces in the neighborhood and return them 
    Halfedge_around_vertex_const_circulator eit = vh->incident_halfedges();  
    do{
      eit++;
      std::copy(eit->surfaces_begin(),eit->surfaces_end(),CGAL::inserter(surfaces)); 
      fh = eit->face(); 
      std::copy(fh->surfaces_begin(),fh->surfaces_end(),CGAL::inserter(surfaces)); 
    }while(eit != vh->incident_halfedges());

    return std::copy(surfaces.begin(),surfaces.end(),oit); 
  }  
  
  CGAL_assertion_msg (false, "Invalid object."); 
  return oit; // never reached  
}


template <typename T>
struct Compare_id 
  : public std::binary_function< T, T, Comparison_result > {
  Comparison_result operator()( const T& x, const T& y) const {
    return CGAL::compare(x.id(),y.id());
  }
};




// returns the first nonnegative root of p
// precondition: p has at least one non negative root
template <class Algebraic_kernel_d_1>
boost::optional< typename Algebraic_kernel_d_1::Algebraic_real_1 >
compute_smallest_nonnegative_root(
    const Algebraic_kernel_d_1& ak,
    const typename Algebraic_kernel_d_1::Polynomial_1& p){
  
  typedef Algebraic_kernel_d_1 AK;
  typedef typename AK::Algebraic_real_1 Root;
  typedef boost::optional< Root > Root_option; 
  
  typename AK::Solve_1 solve_1 = ak.solve_1_object();
  std::vector<Root> roots;
  
  // if p is the zero polynomial then the smallest non negative root is 0 .-)
  if(CGAL::is_zero(p)) return Root_option(Root(0));
  
  solve_1(p,std::back_inserter(roots),false);  

  typename std::vector<Root>::const_iterator it = roots.begin();
  while(it != roots.end() && CGAL::is_negative(*it))
    it++;

  if (it == roots.end())
    return Root_option();
  else 
    return Root_option(*it);
}

template <class Algebraic_kernel_d_1>
CGAL::Comparison_result
compare_smallest_nonnegative_roots(
    const Algebraic_kernel_d_1& ak,
    const typename Algebraic_kernel_d_1::Polynomial_1& p1,
    const typename Algebraic_kernel_d_1::Polynomial_1& p2){
  typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Root; 
  typedef boost::optional<Root> Root_option; 
  Root_option r1  = compute_smallest_nonnegative_root(ak,p1);
  Root_option r2  = compute_smallest_nonnegative_root(ak,p2);
  if(!r1 && !r2 ) return CGAL::EQUAL;
  if(!r1)         return CGAL::LARGER;
  if(!r2)         return CGAL::SMALLER; 
  return                 CGAL::compare(*r1,*r2);
}


} // end of namespace VDOL_3
} //namespace CGAL

#endif // CGAL_VOL_TOOLS_H
