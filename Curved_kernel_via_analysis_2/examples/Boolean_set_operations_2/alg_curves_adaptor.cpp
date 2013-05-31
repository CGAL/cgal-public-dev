/*! \file alg_curves_traits_adapter.cpp
 * Using the traits adaptor to generate a traits class for polygons 
 * bounded by semi-algebraic curves
 */

#include <CGAL/Algebraic_kernel_d/flags.h>

#include <CGAL/basic.h>

#include <iostream>

#ifndef CGAL_USE_CORE

int main ()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl;
  return (0);
}

#else

// Allows to use the Filtered_curve_kernel_via_analysis_2
#ifndef CGAL_ACK_USE_FILTERED_CKvA_2
#define CGAL_ACK_USE_FILTERED_CKvA_2 0
#endif

// What is the coefficient type of the input?
#ifndef CGAL_ACK_COEFFICIENT
#define CGAL_ACK_COEFFICIENT CGAL::CORE_arithmetic_kernel::Integer
#define CGAL_ACK_COEFFICIENT_IS_INTEGER 1
#endif

#include <CGAL/Timer.h>

#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_curve_kernel_2_generator.h>

//#include "include/CGAL/Polynomial_parser_2.h"

#if CGAL_ACK_USE_FILTERED_CKvA_2
#include <CGAL/Filtered_algebraic_curve_kernel_2.h>
#include <CGAL/Curved_kernel_via_analysis_2/Filtered_curved_kernel_via_analysis_2_impl.h>
#else
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h>
#endif

#include <CGAL/Arrangement_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <fstream>

typedef CGAL_ACK_COEFFICIENT Coefficient;

#if !CGAL_ACK_USE_FILTERED_CKvA_2
typedef CGAL::Algebraic_curve_kernel_2_generator<Coefficient>
::Algebraic_curve_kernel_with_qir_and_bitstream_2
Algebraic_curve_kernel_2;
#else
typedef CGAL::Algebraic_curve_kernel_2_generator<Coefficient>
::Filtered_algebraic_curve_kernel_with_qir_and_bitstream_2
Algebraic_curve_kernel_2;
#endif

#if !CGAL_ACK_USE_FILTERED_CKvA_2
typedef CGAL::Curved_kernel_via_analysis_2< Algebraic_curve_kernel_2 > 
Curved_kernel_2; 
#else
typedef CGAL::Curved_kernel_via_analysis_2< Algebraic_curve_kernel_2 > 
Exact_curved_kernel_2; 
typedef CGAL::Filtered_curved_kernel_via_analysis_2<Exact_curved_kernel_2>
Curved_kernel_2; 
#endif

typedef Algebraic_curve_kernel_2::Polynomial_2          Polynomial_2;
typedef Curved_kernel_2                                 Traits_2; 
typedef Traits_2::Curve_2                               Bezier_curve_2;
typedef Traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef Traits_2::Curve_2                               Curve_2;

typedef CGAL::Arrangement_2< Traits_2 >                 Arrangement_2;

typedef CGAL::Gps_traits_2<Traits_2>                    Gps_traits_2;
typedef Gps_traits_2::General_polygon_2                 Polygon_2;
typedef Gps_traits_2::General_polygon_with_holes_2      Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2>                 Polygon_set;

/*! Read a general polygon with holes, formed by semi-algebraic curves, 
 * from the given input file.
 */
bool read_alg_curves_polygon (const char* filename, Polygon_with_holes_2& P)
{
  // Open the input file.
  std::ifstream   in_file (filename);

  if (! in_file.is_open())
    return false;

#if 0
// TODO replace reading
  
  // Read the number of curves.
  unsigned int      n_curves;
  unsigned int      k;

  in_file >> n_curves;

  // Read the curves one by one, and construct the general polygon these
  // curve form (the outer boundary and the holes inside it).
  Traits_2                       traits;
  Traits_2::Make_x_monotone_2    make_x_monotone = 
                                        traits.make_x_monotone_2_object();
  bool                           first = true;
  Rat_point_2                    p_0;
  std::list<X_monotone_curve_2>  xcvs;
  Rat_kernel                     ker;
  Rat_kernel::Equal_2            equal = ker.equal_2_object();
  std::list<Polygon_2>           pgns;

  for (k = 0; k < n_curves; k++) {
    // Read the current curve and subdivide it into x-monotone subcurves.
    Bezier_curve_2                           B;
    std::list<CGAL::Object>                  x_objs;
    std::list<CGAL::Object>::const_iterator  xoit;
    X_monotone_curve_2                       xcv;

    in_file >> B;
    make_x_monotone (B, std::back_inserter (x_objs));
    
    for (xoit = x_objs.begin(); xoit != x_objs.end(); ++xoit) {
      if (CGAL::assign (xcv, *xoit))
        xcvs.push_back (xcv);
    }
    
    // Check if the current curve closes a polygon, namely whether it target
    // point (the last control point) equals the source of the first curve in
    // the current chain.
    if (! first) {
      if (equal (p_0, B.control_point(B.number_of_control_points() - 1))) {
        // Push a new polygon into the polygon list. Make sure that the polygon
        // is counterclockwise oriented if it represents the outer boundary
        // and clockwise oriented if it represents a hole.
        Polygon_2          pgn (xcvs.begin(), xcvs.end());
        CGAL::Orientation  orient = pgn.orientation();
        
        if ((pgns.empty() && orient == CGAL::CLOCKWISE) ||
            (! pgns.empty() && orient == CGAL::COUNTERCLOCKWISE))
          pgn.reverse_orientation();
        
        pgns.push_back (pgn);
        xcvs.clear();
        first = true;
      }
    }
    else {
      // This is the first curve in the chain - store its source point.
      p_0 = B.control_point(0);
      first = false;
    }
  }

  if (! xcvs.empty())
    return false;

  // Construct the polygon with holes.
  std::list<Polygon_2>::iterator     pit = pgns.begin();
  
  ++pit;
  P = Polygon_with_holes_2 (pgns.front(), pit, pgns.end());

#else

  P = Polygon_with_holes_2 ();

#endif
  
  return true;
}


/*! Read a general polygon with holes, formed by semi-algebraic curves, 
 * from the given input file.
 */
bool construct_polygon (int id, Polygon_with_holes_2& P)
{

  Traits_2 traits;
  
  std::list<X_monotone_curve_2> xcvs;
  std::list<Polygon_2> pgns;
  
  Polynomial_2 poly;
  typedef CGAL::Polynomial_traits_d< Polynomial_2 > PT_2;

  PT_2::Construct_polynomial construct_polynomial;
  
  std::list<std::pair<CGAL::Exponent_vector, Coefficient > > icoeffs;

  // define polynomials by id
  if (id == 1) {
    // unit circle
    icoeffs.push_back(std::make_pair(CGAL::Exponent_vector(2,0),1));
    icoeffs.push_back(std::make_pair(CGAL::Exponent_vector(0,2),1));
    icoeffs.push_back(std::make_pair(CGAL::Exponent_vector(0,0),-1));
  } else {
    // axis-aligned ellipse
    icoeffs.push_back(std::make_pair(CGAL::Exponent_vector(2,0),6));
    icoeffs.push_back(std::make_pair(CGAL::Exponent_vector(0,2),2));
    icoeffs.push_back(std::make_pair(CGAL::Exponent_vector(0,0),-4));
  }
  poly  = 
    construct_polynomial(icoeffs.begin(),icoeffs.end());
  std::cout << "The bivariate polynomial: " << poly << std::endl;
  
  Curve_2 curve(poly);
  
  Arrangement_2 arr;
  CGAL::insert(arr, curve);
  
  // iterate faces
  Arrangement_2::Face_iterator fit;
  for (fit = arr.faces_begin(); fit != arr.faces_end(); fit++) {
    // find bounded one
    if (!fit->is_unbounded()) {
      break;
    }
  }
  if (fit == arr.faces_end()) {
    std::cerr << "Algebraic curve with id=" << id << " has no bounded face."
              << std::endl;
    return false;
  }
  // iterate outer ccb
  Arrangement_2::Ccb_halfedge_circulator start, circ;
  circ = start = fit->outer_ccb();
  
  do {

    // set direction
    //std::cout << "Curve: " << circ->curve() << std::flush;
    if (circ->direction() == CGAL::ARR_LEFT_TO_RIGHT) {
      xcvs.push_back(circ->curve());
      //std::cout << " (L2R)"<< std::endl;
    } else {
      xcvs.push_back(traits.construct_opposite_2_object()(circ->curve()));
      //std::cout << " (R2L)"<< std::endl;
    }

    circ++;
    
  } while (circ != start);

  // Push a new polygon into the polygon list. Make sure that the polygon
  // is counterclockwise oriented if it represents the outer boundary
  // and clockwise oriented if it represents a hole.
  Polygon_2 pgn(xcvs.begin(), xcvs.end());
  pgns.push_back(pgn);
  
  // Construct the polygon with holes.
  std::list<Polygon_2>::iterator pit = pgns.begin();
  ++pit;

  P = Polygon_with_holes_2(pgns.front(), pit, pgns.end());
  
  return true;
}


// The main program.
int main (int argc, char **argv)
{

  CGAL::set_pretty_mode(std::cout);

  // Get the name of the input files from the command line, or use the default
  // char_g.dat and char_m.dat files if no command-line parameters are given.
  const char           *filename1 = (argc > 1) ? argv[1] : "poly1.dat";
  const char           *filename2 = (argc > 2) ? argv[2] : "poly2.dat";

  // Read the general polygons from the input files.
  CGAL::Timer           timer;

  Polygon_with_holes_2  P1, P2;

  timer.start();

  if (! read_alg_curves_polygon (filename1, P1)) {
    std::cerr << "Failed to read P1: " << filename1 << " ..." << std::endl;
    return 1;
  }

  if (! read_alg_curves_polygon (filename2, P2)) {
    std::cerr << "Failed to read P2: " << filename2 << " ..." << std::endl;
    return 1;
  }

  // read a bounded face (construct arrangement)
  construct_polygon(1, P1);
  construct_polygon(2, P2);
  
  timer.stop();
  std::cout << "Constructed the input polygons in " << timer.time() 
            << " seconds." << std::endl << std::endl;
  
  // Compute the intersection of the two polygons.
  Polygon_set                     I, S, J;
  Polygon_set::const_iterator     psit;

  /////////////////////////////////////////////////////////////////////////////
  timer.reset();
  timer.start();
  CGAL::intersection (P1, P2, std::back_inserter(I));
  timer.stop();

  std::cout << "The intersection polygons are of sizes: {";
  for (psit = I.begin(); psit != I.end(); ++psit) {
    std::cout << " " << psit->outer_boundary().size() << " (h=";
    std::cout << psit->number_of_holes() << ")";
  }
  std::cout << " }" << std::endl;
  std::cout << "The intersection computation took "
            << timer.time() << " seconds." << std::endl;
  std::cout << std::endl;

  /////////////////////////////////////////////////////////////////////////////
  timer.reset();
  timer.start();
  CGAL::symmetric_difference (P1, P2, std::back_inserter(S));
  timer.stop();

  std::cout << "The symmetric difference polygons are of sizes: {";
  for (psit = S.begin(); psit != S.end(); ++psit) {
    std::cout << " " << psit->outer_boundary().size() << " (h=";
    std::cout << psit->number_of_holes() << ")";
  }
  std::cout << " }" << std::endl;
  std::cout << "The symmetric difference computation took "
            << timer.time() << " seconds." << std::endl;
  std::cout << std::endl;

  /////////////////////////////////////////////////////////////////////////////
#if 1
  // TODO join requires aggregrated construction - BSO_2 does currently
  //      not supported aggregated constructions for "open" topologies
  //      More detailed: reference_face is not yet implemented
  timer.reset();
  timer.start();
  std::list< Polygon_with_holes_2 > polygons;
  polygons.push_back(P1);
  polygons.push_back(P2);
  CGAL::join (polygons.begin(), polygons.end(), std::back_inserter(J));
  timer.stop();

  std::cout << "The join polygons are of sizes: {";
  for (psit = J.begin(); psit != J.end(); ++psit) {
    std::cout << " " << psit->outer_boundary().size() << " (H=";
    std::cout << psit->number_of_holes() << ")";
  }
  std::cout << " }" << std::endl;
  std::cout << "The join computation took "
            << timer.time() << " seconds." << std::endl;
  std::cout << std::endl;
#endif

  return 0;
}

#endif
