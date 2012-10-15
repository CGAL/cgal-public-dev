// this is to enable the arc ends at infinity in CKvA
#define CGAL_ARRANGEMENT_ON_DUPIN_CYCLIDE 1
#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1 

#define CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE 1
#define CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE 1
#define CGAL_SL_VERBOSE 1

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/macros.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d_1.h>

// consider this once things are running again... 
//#include <CGAL/Cached_algebraic_kernel_1.h>
//#include <CGAL/Cached_algebraic_kernel_2.h>

#include <CGAL/Algebraic_kernel_d/Algebraic_curve_kernel_2.h> // template arg to CKvA_2
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h> // traits for Arr_2
#include <CGAL/Algebraic_curve_kernel_2_generator.h> // generates a default kernel 
#include <CGAL/Envelope_3/Envelope_diagram_on_surface_2.h>

#include <CGAL/VDOL_3/Single_voronoi_cell_envelope_traits_3.h>
#include <CGAL/VDOL_3/Cylinder_traits_3.h>
#include <CGAL/Arr_cylindrical_topology_traits_2.h>
#include <CGAL/envelope_3.h>


#include <CGAL/construct_point_2.h>


// get default arithmetic 
typedef CGAL::Arithmetic_kernel AK;
typedef AK::Integer Integer;
typedef AK::Rational Rational;

// define linear kernel 
typedef CGAL::Cartesian<Rational>                       Linear_kernel_3;
typedef Linear_kernel_3::Line_3                         Line_3;
typedef Linear_kernel_3::Point_3                        Point_3;

// define algebraic kerenel 
typedef CGAL::Algebraic_kernel_d_1_generator<Integer>
::Default_algebraic_kernel_1                            AK_1;

typedef CGAL::Algebraic_curve_kernel_2<AK_1>            AK_2;

// define curved kernel 
typedef CGAL::Curved_kernel_via_analysis_2< AK_2 >      Curved_kernel_2;

// define traits for lower envelop diagram 
typedef CGAL::VDOL_3::Single_voronoi_cell_envelope_traits_3
<Linear_kernel_3, Curved_kernel_2>                 SVCET_3;


typedef CGAL::VDOL_3::Cylinder_traits_3<SVCET_3> CT_3; // Geometry_traits 
typedef CGAL::Envelope_3::Envelope_pm_dcel<CT_3,CT_3::Xy_monotone_surface_3> Envelope_pm_dcel;

typedef CGAL::Arr_cylindrical_topology_traits_2<CT_3,Envelope_pm_dcel> CTT_3;

typedef CGAL::Envelope_diagram_on_surface_2<CT_3,CTT_3> Envelop_diagram_on_cylinder_3;




int main(int argc, char* argv[])
{ 

  Line_3 line_0 (Point_3( 0, 0, 0), Point_3( 1, 0, 0));
  Line_3 line_1(Point_3(0,0,1),Point_3(0,1,1));
  
  std::list<Line_3> lines;
  lines.push_back(line_1);

  boost::shared_ptr<CT_3> ct_3 =
    boost::make_shared<CT_3>(CT_3(line_0,4));
  
  Envelop_diagram_on_cylinder_3 diagram_on_cylinder(&(*ct_3));
  CGAL::lower_envelope_3(lines.begin(), lines.end(), diagram_on_cylinder);

  const SVCET_3* svcet_3 = ct_3->svcet_3();
  CT_3::Point_2 nbpoint_1(construct_point_2(svcet_3,Rational(0),CGAL::ARR_BOTTOM_BOUNDARY),CGAL::NEGATIVE);
  CT_3::Point_2 nbpoint_2(construct_point_2(svcet_3,Rational(1),CGAL::ARR_BOTTOM_BOUNDARY),CGAL::NEGATIVE);

  CT_3::Point_2 nipoint_1(construct_point_2(svcet_3,Rational(0),Rational(0)),CGAL::NEGATIVE);
  CT_3::Point_2 nipoint_2(construct_point_2(svcet_3,Rational(1),Rational(0)),CGAL::NEGATIVE);
  CT_3::Point_2 nipoint_3(construct_point_2(svcet_3,Rational(1),Rational(1)),CGAL::NEGATIVE);

  CT_3::Point_2 ntpoint_1(construct_point_2(svcet_3,Rational(0),CGAL::ARR_TOP_BOUNDARY),CGAL::NEGATIVE);
  CT_3::Point_2 ntpoint_2(construct_point_2(svcet_3,Rational(1),CGAL::ARR_TOP_BOUNDARY),CGAL::NEGATIVE);

  CT_3::Point_2 pbpoint_1(construct_point_2(svcet_3,Rational(0),CGAL::ARR_BOTTOM_BOUNDARY),CGAL::POSITIVE);
  CT_3::Point_2 pbpoint_2(construct_point_2(svcet_3,Rational(1),CGAL::ARR_BOTTOM_BOUNDARY),CGAL::POSITIVE);

  CT_3::Point_2 pipoint_1(construct_point_2(svcet_3,Rational(0),Rational(0)),CGAL::POSITIVE);
  CT_3::Point_2 pipoint_2(construct_point_2(svcet_3,Rational(1),Rational(0)),CGAL::POSITIVE);
  CT_3::Point_2 pipoint_3(construct_point_2(svcet_3,Rational(1),Rational(1)),CGAL::POSITIVE);

  CT_3::Point_2 ptpoint_1(construct_point_2(svcet_3,Rational(0),CGAL::ARR_TOP_BOUNDARY),CGAL::POSITIVE);
  CT_3::Point_2 ptpoint_2(construct_point_2(svcet_3,Rational(1),CGAL::ARR_TOP_BOUNDARY),CGAL::POSITIVE);
  
  {
    //TEST Compare_x_2 
    // this is just a forwards since planes are on top of each other
    CT_3::Compare_x_2 compare_x_2 = ct_3->compare_x_2_object();
    assert(compare_x_2(nbpoint_1,nbpoint_1) == CGAL::EQUAL);
    assert(compare_x_2(nbpoint_1,nbpoint_2) == CGAL::SMALLER);
    assert(compare_x_2(nbpoint_1,nipoint_1) == CGAL::EQUAL);
    assert(compare_x_2(nbpoint_1,nipoint_2) == CGAL::SMALLER);
    assert(compare_x_2(nbpoint_1,nipoint_3) == CGAL::SMALLER);
    assert(compare_x_2(nbpoint_1,ntpoint_1) == CGAL::EQUAL);
    assert(compare_x_2(nbpoint_1,ntpoint_2) == CGAL::SMALLER);
    assert(compare_x_2(nbpoint_1,pbpoint_1) == CGAL::EQUAL);
    assert(compare_x_2(nbpoint_1,pbpoint_2) == CGAL::SMALLER);
    assert(compare_x_2(nbpoint_1,pipoint_1) == CGAL::EQUAL);
    assert(compare_x_2(nbpoint_1,pipoint_2) == CGAL::SMALLER);
    assert(compare_x_2(nbpoint_1,pipoint_3) == CGAL::SMALLER);
    assert(compare_x_2(nbpoint_1,ptpoint_1) == CGAL::EQUAL);
    assert(compare_x_2(nbpoint_1,ptpoint_2) == CGAL::SMALLER);
    
    assert(compare_x_2(pipoint_2,nbpoint_1) == CGAL::LARGER);
    assert(compare_x_2(pipoint_2,nbpoint_2) == CGAL::EQUAL);
    assert(compare_x_2(pipoint_2,nipoint_1) == CGAL::LARGER);
    assert(compare_x_2(pipoint_2,nipoint_2) == CGAL::EQUAL);
    assert(compare_x_2(pipoint_2,nipoint_3) == CGAL::EQUAL);
    assert(compare_x_2(pipoint_2,ntpoint_1) == CGAL::LARGER);
    assert(compare_x_2(pipoint_2,ntpoint_2) == CGAL::EQUAL);
    assert(compare_x_2(pipoint_2,pbpoint_1) == CGAL::LARGER);
    assert(compare_x_2(pipoint_2,pbpoint_2) == CGAL::EQUAL);
    assert(compare_x_2(pipoint_2,pipoint_1) == CGAL::LARGER);
    assert(compare_x_2(pipoint_2,pipoint_2) == CGAL::EQUAL);
    assert(compare_x_2(pipoint_2,pipoint_3) == CGAL::EQUAL);
    assert(compare_x_2(pipoint_2,ptpoint_1) == CGAL::LARGER);
    assert(compare_x_2(pipoint_2,ptpoint_2) == CGAL::EQUAL);
  }
  {    
    CT_3::Compare_xy_2 compare_xy_2 = ct_3->compare_xy_2_object();
    assert(compare_xy_2(nbpoint_1,nipoint_1) == CGAL::LARGER);
    assert(compare_xy_2(nbpoint_1,ntpoint_1) == CGAL::LARGER);
    assert(compare_xy_2(nbpoint_1,pbpoint_1) == CGAL::LARGER);
    assert(compare_xy_2(nbpoint_1,pipoint_1) == CGAL::LARGER);
    assert(compare_xy_2(nbpoint_1,ptpoint_1) == CGAL::EQUAL);
    assert(compare_xy_2(nbpoint_1,nbpoint_1) == CGAL::EQUAL);
    
    assert(compare_xy_2(nbpoint_1,nipoint_2) == CGAL::SMALLER);
    assert(compare_xy_2(nbpoint_1,ntpoint_2) == CGAL::SMALLER);
    assert(compare_xy_2(nbpoint_1,pbpoint_2) == CGAL::SMALLER);
    assert(compare_xy_2(nbpoint_1,pipoint_2) == CGAL::SMALLER);
    assert(compare_xy_2(nbpoint_1,ptpoint_2) == CGAL::SMALLER);
    assert(compare_xy_2(nbpoint_1,nbpoint_2) == CGAL::SMALLER);
    
    assert(compare_xy_2(nbpoint_1,nipoint_3) == CGAL::SMALLER);
    assert(compare_xy_2(nbpoint_1,pipoint_3) == CGAL::SMALLER);
   
    assert(compare_xy_2(pipoint_2,nipoint_1) == CGAL::LARGER);
    assert(compare_xy_2(pipoint_2,ntpoint_1) == CGAL::LARGER);
    assert(compare_xy_2(pipoint_2,pbpoint_1) == CGAL::LARGER);
    assert(compare_xy_2(pipoint_2,pipoint_1) == CGAL::LARGER);
    assert(compare_xy_2(pipoint_2,ptpoint_1) == CGAL::LARGER);
    assert(compare_xy_2(pipoint_2,nbpoint_1) == CGAL::LARGER);
    
    assert(compare_xy_2(pipoint_2,nipoint_2) == CGAL::LARGER);
    assert(compare_xy_2(pipoint_2,ntpoint_2) == CGAL::LARGER);
    assert(compare_xy_2(pipoint_2,pbpoint_2) == CGAL::LARGER);
    assert(compare_xy_2(pipoint_2,pipoint_2) == CGAL::EQUAL);
    assert(compare_xy_2(pipoint_2,ptpoint_2) == CGAL::SMALLER);
    assert(compare_xy_2(pipoint_2,nbpoint_2) == CGAL::SMALLER);
    
    assert(compare_xy_2(pipoint_2,nipoint_3) == CGAL::LARGER);
    assert(compare_xy_2(pipoint_2,pipoint_3) == CGAL::SMALLER);


    assert(compare_xy_2(nbpoint_1,nipoint_1) == CGAL::LARGER);
    assert(compare_xy_2(nbpoint_1,ntpoint_1) == CGAL::LARGER);
    assert(compare_xy_2(nbpoint_1,pbpoint_1) == CGAL::LARGER);
    assert(compare_xy_2(nbpoint_1,pipoint_1) == CGAL::LARGER);
    assert(compare_xy_2(nbpoint_1,ptpoint_1) == CGAL::EQUAL);
    assert(compare_xy_2(nbpoint_1,nbpoint_1) == CGAL::EQUAL);

    assert(compare_xy_2(nipoint_1,nipoint_1) == CGAL::EQUAL);
    assert(compare_xy_2(nipoint_1,ntpoint_1) == CGAL::SMALLER);
    assert(compare_xy_2(nipoint_1,pbpoint_1) == CGAL::SMALLER);
    assert(compare_xy_2(nipoint_1,pipoint_1) == CGAL::SMALLER);
    assert(compare_xy_2(nipoint_1,ptpoint_1) == CGAL::SMALLER);
    assert(compare_xy_2(nipoint_1,nbpoint_1) == CGAL::SMALLER);

    assert(compare_xy_2(ntpoint_1,nipoint_1) == CGAL::LARGER);
    assert(compare_xy_2(ntpoint_1,ntpoint_1) == CGAL::EQUAL);
    assert(compare_xy_2(ntpoint_1,pbpoint_1) == CGAL::EQUAL);
    assert(compare_xy_2(ntpoint_1,pipoint_1) == CGAL::SMALLER);
    assert(compare_xy_2(ntpoint_1,ptpoint_1) == CGAL::SMALLER);
    assert(compare_xy_2(ntpoint_1,nbpoint_1) == CGAL::SMALLER);

    assert(compare_xy_2(ntpoint_1,nipoint_1) == CGAL::LARGER);
    assert(compare_xy_2(ntpoint_1,ntpoint_1) == CGAL::EQUAL);    
    assert(compare_xy_2(ntpoint_1,pbpoint_1) == CGAL::EQUAL);
    assert(compare_xy_2(ntpoint_1,pipoint_1) == CGAL::SMALLER);
    assert(compare_xy_2(ntpoint_1,ptpoint_1) == CGAL::SMALLER);
    assert(compare_xy_2(ntpoint_1,nbpoint_1) == CGAL::SMALLER);
  }
}

