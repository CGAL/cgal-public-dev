#include <CGAL/basic.h> //Support for standard header names
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/CORE_BigInt.h>
#include <CGAL/Algebraic_kernel_d_1.h> //Algebraic  Kernel, needed by rational Traits.
#include <CGAL/Arrangement_2.h> //Arrangement class
#include <CGAL/Arr_rational_function_traits_2.h> //Rational function traits
#include <CGAL/Arr_naive_point_location.h> // Naive point location to verify results
#include <CGAL/Arr_landmarks_point_location.h> // Main PL strategy to test, now supporting Rational functions
#include <CGAL/Arr_point_location/Arr_lm_random_generator.h> //Landmark Generators
#include <CGAL/Arr_point_location/Arr_lm_grid_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_halton_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_middle_edges_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_specified_points_generator.h>

#include <CGAL/Timer.h> // Timer for benchmarks
#include <CGAL/CORE_BigInt.h>                      // NT
#include <CGAL/Algebraic_kernel_d_1.h>             // Algebraic Kernel
#include <CGAL/Arr_rational_function_traits_2.h>   // Traits
#include <CGAL/Arrangement_2.h>                    // Arrangement
typedef CORE::BigInt Number_type;
typedef CGAL::Algebraic_kernel_d_1<Number_type> AK1;
typedef CGAL::Arr_rational_function_traits_2<AK1> Traits_2;

typedef Traits_2::Polynomial_1 Polynomial_1;
typedef Traits_2::Algebraic_real_1 Alg_real_1;

typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
typedef Arrangement_2::Halfedge_handle Halfedge_handle;
typedef Arrangement_2::Edge_const_iterator Edge_const_iterator;
typedef Arrangement_2::Vertex_const_iterator Vertex_const_iterator;

typedef CGAL::Arr_naive_point_location<Arrangement_2> Naive_point_location;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2> Lm_point_location;
typedef CGAL::Arr_random_landmarks_generator<Arrangement_2> Random_lm_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Random_lm_generator> Lm_random_point_location;
typedef CGAL::Arr_grid_landmarks_generator<Arrangement_2> Grid_lm_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Grid_lm_generator> Lm_grid_point_location;
typedef CGAL::Arr_halton_landmarks_generator<Arrangement_2> Halton_lm_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Halton_lm_generator> Lm_halton_point_location;
typedef CGAL::Arr_middle_edges_landmarks_generator<Arrangement_2> Middle_edges_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Middle_edges_generator> Lm_middle_edges_point_location;
typedef CGAL::Arr_landmarks_specified_points_generator<Arrangement_2> Specified_points_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2,
Specified_points_generator> Lm_specified_points_point_location;

typedef Traits_2::Polynomial_1 Polynomial_1;
typedef Traits_2::Algebraic_real_1 Alg_real_1;
typedef Traits_2::Point_2 Point_2;
typedef Traits_2::Curve_2 Curve_2;

typedef std::list<Point_2>                                Points_list;
typedef Points_list::iterator                             Point_iterator;
typedef std::list<Curve_2>                                Curve_list;
typedef std::vector<CGAL::Object>                         Objects_vector;
typedef Objects_vector::iterator                          Object_iterator;

// ===> Change the number of point-location startegies
//      when a new point location is added    <===
#define NUM_OF_POINT_LOCATION_STRATEGIES 3

/*!
  Prints all x-monotone curves along a given CCB.
  Used for debugging
*/
void print_ccb (Arrangement_2::Ccb_halfedge_const_circulator circ)
{
    Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
    //std::cout << "(" << curr->source()->point() << ")";
    do {

        if(curr->is_fictitious()){
            std::cout << "Fictitious Edge" << std::endl;
            std::cout << curr->curve() << std:: endl;
            continue;
        }

        //Arrangement_2::Halfedge_const_handle he = curr ;
        std::cout << "   [" << curr->curve() << "]   "
                  << "(" << curr->target()->point() << ")";
    } while (++curr != circ);
    std::cout << std::endl;
}

/*!
  Prints all x-monotone curves along a given unbounded CCB.
  Used for debugging.
*/
void print_unbounded_face(Arrangement_2::Face_const_handle f)
{

    Arrangement_2::Ccb_halfedge_const_circulator curr, first;
    Arrangement_2::Halfedge_const_handle he = curr;

    curr = first = f->outer_ccb();
    if (! curr->source()->is_at_open_boundary())
        std::cout << "(" << curr->source()->point() << ")";
    do {
        he = curr;
        if (! he->is_fictitious())
            std::cout << "   [" << he->curve() << "]   ";
        else
            std::cout << "   [ ... ]   ";

        if (! he->target()->is_at_open_boundary())
            std::cout << "(" << he->target()->point() << ")";

        ++curr;
    } while (curr != first);
    std::cout << std::endl;
}

/*!
    Prints the outer and inner boundaries of a given face
*/
void print_face (Arrangement_2::Face_const_handle f)
{
    // Print the outer boundary.
    if (f->is_unbounded()){
        std::cout << "Unbounded face. " << std::endl;
        print_unbounded_face(f);
    }
    else {
        std::cout << "Outer boundary: ";
        print_ccb (f->outer_ccb());
    }

    // Print the boundary of each of the holes.
    Arrangement_2::Hole_const_iterator hi;
    int                                 index = 1;
    for (hi = f->holes_begin(); hi != f->holes_end(); ++hi, ++index) {
        std::cout << "    Hole #" << index << ": ";
        print_ccb (*hi);
    }

    // Print the isolated vertices.
    Arrangement_2::Isolated_vertex_const_iterator iv;
    for (iv = f->isolated_vertices_begin(), index = 1;
         iv != f->isolated_vertices_end(); ++iv, ++index)
    {
        std::cout << "    Isolated vertex #" << index << ": "
                  << "(" << iv->point() << ")" << std::endl;
    }
}


/*!
    Execute PL point location of the given Points_list
    Append resulting objects to Out container
    returns location run time.
*/
template <class PL, class Out>
double run_pl (Points_list &plist, PL &pl, Out out)
{


    Arrangement_2::Vertex_const_handle v;
    Arrangement_2::Halfedge_const_handle e;
    Arrangement_2::Face_const_handle f;

    CGAL::Timer timer;
    timer.reset();
    timer.start(); //START
    for (Point_iterator piter = plist.begin(); piter != plist.end(); piter++) {
        CGAL::Object obj = pl.locate (*piter);
        CGAL_assertion(obj.is_empty() == false);
        *out++ = obj;

        std::cout << "The point " << *piter << " is located ";
        if (CGAL::assign(f, obj))
        {
            // q is located inside a face:
            if (f->is_unbounded())
                std::cout << "inside an unbounded face." << std::endl;
            else
                std::cout << "inside a bounded face." << std::endl;

            print_face(f);
        }
        else if (CGAL::assign(e, obj))
        {
            // q is located on an edge:
            std::cout << "on an edge: " << e->curve() << std::endl;
        }
        else if (CGAL::assign(v, obj))
        {
            // q is located on a vertex:
            if (v->is_isolated())
                std::cout << "on an isolated vertex: " << v->point() << std::endl;
            else
                std::cout << "on a vertex: " << v->point() << std::endl;
        }
        else
        {
            CGAL_assertion_msg(false, "Invalid object.");
        }

    }
    timer.stop(); ///END


    return timer.time();
}

/*!
   Check point location variants (different generators) againts the given Arrangement_2
   using Points_list as query points
 */
int check_point_location (Arrangement_2 &arr, Points_list &plist)
{
    //init - all point locations
    CGAL::Timer            timer;

    //Used to validate results
    Naive_point_location            naive_pl (arr);                 // 0

    //Only LM gen that works at the moment with Rational function Traits
    timer.reset(); timer.start();
    Lm_point_location               lm_pl (arr);                    // 1
    timer.stop();
    std::cout << "Lm (vert) construction took " << timer.time() <<std::endl;

    Specified_points_generator::Points_set points;
    AK1 ak1;
    Traits_2 traits(&ak1);
    Traits_2::Construct_point_2 construct_point = traits.construct_point_2_object();

    //Single point to force walks
    points.push_back(construct_point(CORE::BigInt(-4), CORE::BigInt(-4)));

    timer.reset(); timer.start();
    Specified_points_generator                specified_points_g(arr,points); //2
    Lm_specified_points_point_location        specified_points_lm_pl (arr, &specified_points_g);
    timer.stop();
    std::cout << "Specified_points lm construction took " << timer.time() <<std::endl;

    /*
     timer.reset(); timer.start();
     Random_lm_generator             random_g(arr);
     Lm_random_point_location        random_lm_pl (arr, &random_g);  // 5
     Lm_random_point_location_no_c   random_lm_noc_pl (arr, &random_g);  // 6
     timer.stop();
     std::cout << "Random lm construction took " << timer.time() <<std::endl;

     timer.reset(); timer.start();
     Grid_lm_generator               grid_g(arr);
     Lm_grid_point_location          grid_lm_pl (arr, &grid_g);      // 7
     Lm_grid_point_location_no_c     grid_lm_noc_pl (arr, &grid_g);      // 8
     timer.stop();

     td::cout << "Grid lm construction took " << timer.time() <<std::endl;


     timer.reset(); timer.start();
     Halton_lm_generator             halton_g(arr);
     Lm_halton_point_location        halton_lm_pl (arr, &halton_g);  // 9
     Lm_halton_point_location_no_c   halton_lm_noc_pl (arr, &halton_g);  // 10
     timer.stop();
     std::cout << "Halton lm construction took " << timer.time() <<std::endl;
          timer.reset(); timer.start();
     Middle_edges_generator             middle_edges_g(arr);
     Lm_middle_edges_point_location        middle_edges_lm_pl (arr, &middle_edges_g);  // 11
     Lm_middle_edges_point_location_no_c   middle_edges_lm_noc_pl (arr, &middle_edges_g);  // 12
     timer.stop();
     std::cout << "Middle edges lm construction took " << timer.time() <<std::endl;

     */

    // ===> Add new point location instance here. <===

    Objects_vector                  objs[NUM_OF_POINT_LOCATION_STRATEGIES];
    Object_iterator                 ob_iter[NUM_OF_POINT_LOCATION_STRATEGIES];
    Arrangement_2::Vertex_const_handle    vh_ref, vh_curr;
    Arrangement_2::Halfedge_const_handle  hh_ref, hh_curr;
    Arrangement_2::Face_const_handle      fh_ref, fh_curr;

    //LOCATE the points in the list using all PL strategies

    //std::cout << "Time in seconds" <<std::endl; ;
    std::cout << std::endl;

    std::cout << "Naive location took " << run_pl(plist, naive_pl, std::back_inserter(objs[0])) << std::endl;
    std::cout << "Landmarks (vertices) location took "  << run_pl(plist, lm_pl, std::back_inserter(objs[1])) << std::endl;
    std::cout << "Landmarks (specified) location took "  << run_pl(plist, specified_points_lm_pl, std::back_inserter(objs[2])) << std::endl;

    //END LOCATION
    int pls_num = NUM_OF_POINT_LOCATION_STRATEGIES;
    int pl_index;
    int result = 0;

    //Init all obejct iterators
    for (pl_index=0; pl_index<pls_num; pl_index++)
    {
        ob_iter[pl_index] = objs[pl_index].begin();
    }

    //get size of objects
    unsigned int size = objs[0].size();
    //std::cout <<"size is "<< size << std::endl;

    for (pl_index=0; pl_index<pls_num; pl_index++)
    {
        if (size != objs[pl_index].size())
        {
            std::cout << "Error: size of pl number "<<pl_index<<" is "
                      <<objs[pl_index].size()<< std::endl;
            result = -1;
        }
    }

    //Assign and check results are the same in all PL

    unsigned int qi; //qi is the query point index

    for (qi=0; qi<size; qi++)
    {
        //assign object to a face
        if (CGAL::assign (fh_ref, ob_iter[0][qi]))
        {
            for (int pl_index=1; pl_index<pls_num; pl_index++)
            {
                if (! CGAL::assign(fh_curr, ob_iter[pl_index][qi]))
                {
                    std::cout << "Error in point location number " << pl_index;
                    if (CGAL::assign(fh_curr, ob_iter[pl_index][qi]))
                    {
                        std::cout << ", an halfedge returned instead of a face"<< std::endl;
                    }
                    else if (CGAL::assign(hh_curr, ob_iter[pl_index][qi]))
                    {
                        std::cout << ", a vertex returned instead of a face"<< std::endl;
                    }
                    else
                    {
                        std::cout << ", an unknown object returned instead of a face"<< std::endl;
                    }
                    result = -1;
                }
                else if (fh_curr != fh_ref)
                {
                    std::cout << "Error: point location number "
                              << pl_index << " return a different face"<< std::endl;
                    result = -1;
                }
            }
        }

        //assign object to a halfedge
        else if (CGAL::assign (hh_ref, ob_iter[0][qi]))
        {
            std::cout << "Halfedge: "<< hh_ref->curve() << std::endl;
            for (int pl_index=1; pl_index<pls_num; pl_index++)
            {
                if (! CGAL::assign(hh_curr, ob_iter[pl_index][qi]))
                {
                    std::cout << "Error in point location number " << pl_index;
                    if (CGAL::assign(fh_curr, ob_iter[pl_index][qi]))
                    {
                        std::cout << ", a face returned instead of an halfedge"<< std::endl;
                    }
                    else if (CGAL::assign(hh_curr, ob_iter[pl_index][qi]))
                    {
                        std::cout << ", a vertex returned instead of an halfedge"<< std::endl;
                    }
                    else
                    {
                        std::cout << ", an unknown object returned instead of an halfedge"<< std::endl;
                    }
                    result = -1;
                }
                else if ((hh_curr != hh_ref) && (hh_curr->twin() != hh_ref))
                {
                    std::cout << "Error: point location number "
                              << pl_index << " return a different halfedge"<< std::endl;
                    std::cout << "Halfedge (curr): "<< hh_curr->curve() << std::endl;
                    result = -1;
                }
            }
        }

        //assign object to a vertex
        else if (CGAL::assign (vh_ref, ob_iter[0][qi]))
        {
            for (int pl_index=1; pl_index<pls_num; pl_index++)
            {
                if (! CGAL::assign(vh_curr, ob_iter[pl_index][qi]))
                {
                    std::cout << "Error in point location number " << pl_index;
                    if (CGAL::assign(fh_curr, ob_iter[pl_index][qi]))
                    {
                        std::cout << ", a face returned instead of a vertex"<< std::endl;
                    }
                    else if (CGAL::assign(hh_curr, ob_iter[pl_index][qi]))
                    {
                        std::cout << ", an halfedge returned instead of a vertex"<< std::endl;
                    }
                    else
                    {
                        std::cout << ", an unknown object returned instead of a vertex"<< std::endl;
                    }
                    result = -1;
                }
                else if (vh_curr != vh_ref)
                {
                    std::cout << "Error: point location number "
                              << pl_index << " return a different vertex"<< std::endl;
                    result = -1;
                }
            }
            std::cout << "Vertex: "<< vh_ref->point() << std::endl;
        }

        else
        {
            std::cout << "Illegal point-location result." << std::endl;
            result = -1;
        }
    }

    return (result);
}

/*!
    Inserts curves to the arrangement and measures insertion time
*/
void insert_curves(Arrangement_2 &arr, Curve_list &curve_list)
{

    CGAL::Timer timer;
    timer.reset();
    timer.start(); //START
    insert (arr, curve_list.begin(), curve_list.end());
    timer.stop(); ///END
    std::cout << "Arrangement aggregate construction took "
              << timer.time() <<std::endl;
    // Print the size of the arrangement.
    std::cout << "V = " << arr.number_of_vertices()
              << ",  E = " << arr.number_of_edges()
              << ",  F = " << arr.number_of_faces() << std::endl;
}

//-TESTS---------------------------------------------------------------------------

AK1 ak1;
Traits_2 traits(&ak1);
Traits_2::Construct_curve_2 construct = traits.construct_curve_2_object();
// a polynomial representing x .-)
Polynomial_1 x = CGAL::shift(Polynomial_1(1), 1);
// 1 and -1
Polynomial_1 P1(1);
Polynomial_1 minusP1(-P1);

/*!
 Empty arrangement
*/
bool test0(){

    //-Initialization-------------------------------------------------
    Arrangement_2 arr;
    Curve_list curve_list;
    Points_list plist;
    //-Curves definition-----------------------------------------------------------
    // no curves
    insert_curves(arr, curve_list);
    //-Point generation------------------------------------------------------------
    plist.push_back(traits.construct_point_2_object()(CORE::BigRat(0), CORE::BigInt(0)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(0), CORE::BigInt(3)));
    //-----------------------------------------------------------------------------

    //check point location of points
    if (check_point_location(arr, plist))
    {
        std::cout << "ERROR in check_point_location."<<std::endl<<std::endl;
        return (false);
    }
    std::cout << std::endl;
    return (true);
}

/*!
   Example from User Manual
*/
bool test1()
{

    //-Initialization-------------------------------------------------
    Arrangement_2 arr;
    Curve_list curve_list;
    Points_list plist;

    //-Curves definition----------------------------------------------
    // container storing all arcs
    std::vector<Traits_2::Curve_2> arcs;

    // Create the rational functions (y = 1 / x), and (y = -1 / x).
    Polynomial_1 Q1 = x;
    Polynomial_1 P0(0);

    Curve_2 c_P1(construct(P1, Q1));
    curve_list.push_back(c_P1);
    curve_list.push_back(construct(minusP1, Q1));

    // Create a bounded segments of the parabolas (y = -4*x^2 + 3) and
    // (y = 4*x^2 - 3), defined over [-sqrt(3)/2, sqrt(3)/2].
    Polynomial_1 P2 = -4 * x * x + 3;
    Polynomial_1 minusP2 = -P2;
    std::vector<std::pair<Alg_real_1, int> > roots;

    // [-sqrt(3)/2, sqrt(3)/2]
    traits.algebraic_kernel_d_1()->solve_1_object()(P2,
                                                    std::back_inserter(roots));
    curve_list.push_back(construct(P2, roots[0].first, roots[1].first));
    curve_list.push_back(construct(minusP2, roots[0].first, roots[1].first));

    // Create the rational function (y = 1 / 2*x) for x > 0, and the
    // rational function (y = -1 / 2*x) for x < 0.
    Polynomial_1 P3(1);
    Polynomial_1 minusP3(-P3);
    Polynomial_1 Q3 = 2 * x;

    curve_list.push_back(construct(P3, Q3, Alg_real_1(0), true));
    curve_list.push_back(construct(minusP3, Q3, Alg_real_1(0), false));

    insert_curves(arr, curve_list);

    //-Point generation------------------------------------------------------------
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(0), CORE::BigInt(0)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(0), CORE::BigInt(3)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigRat(0.17), CORE::BigRat(2.88)));
    //-----------------------------------------------------------------------------

    //check point location of points
    if (check_point_location(arr, plist))
    {
        std::cout << "ERROR in check_point_location."<<std::endl<<std::endl;
        return (false);
    }
    std::cout << std::endl;
    return (true);
}


/*!

*/
bool test2()
{

    //-Initialization-------------------------------------------------
    Arrangement_2 arr;
    Curve_list curve_list;
    Points_list plist;
    //-Curves definition----------------------------------------------
    /*
     y = x^2-4,
     y = -x^2+4
    */
    Polynomial_1 x2_4 = (x*x)-4;
    Polynomial_1 x2plus4 = -(x*x)+4;

    curve_list.push_back(construct(x2_4,P1));
    curve_list.push_back(construct(x2plus4,P1));

    insert_curves(arr, curve_list);

    //-Point generation------------------------------------------------------------
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(0), CORE::BigInt(0)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(2), CORE::BigInt(0)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(-2), CORE::BigInt(0)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(0), CORE::BigInt(4)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(0), CORE::BigInt(-4)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(1), CORE::BigInt(-3)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(1), CORE::BigInt(5)));
    //-----------------------------------------------------------------------------

    //check point location of points
    if (check_point_location(arr, plist))
    {
        std::cout << "ERROR in check_point_location."<<std::endl<<std::endl;
        return (false);
    }
    std::cout << std::endl;

    return (true);
}

/*!
  */
bool test3()
{

    //-Initialization-------------------------------------------------
    Arrangement_2 arr;
    Curve_list curve_list;
    Points_list plist;
    //-Curves definition----------------------------------------------
    /*
     xy=1,
     xy=2,
     xy=3,
     y = (x-4)^2+3,
     y = -(x-4)^2 + 4
    */

    // xy=1
    curve_list.push_back(construct(P1,x));
    // xy=2
    curve_list.push_back(construct(P1*2,x));
    // xy=3
    curve_list.push_back(construct(P1*3,x));
    // y = (x-4)^2+3 [3,5]
    curve_list.push_back(construct((x-4)*(x-4)+3,P1, Alg_real_1(3), Alg_real_1(5)));
    //y = -(x-4)^2 + 4 [3,5]
    curve_list.push_back(construct(minusP1*(x-4)*(x-4)+4,P1, Alg_real_1(3), Alg_real_1(5)));

    insert_curves(arr, curve_list);

    //-Point generation------------------------------------------------------------
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(0), CORE::BigInt(0)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(5), CORE::BigInt(5)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigRat(1.3), CORE::BigRat(1)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigRat(1.5), CORE::BigRat(1.5)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigRat(3.5), CORE::BigRat(3.5)));
    //-----------------------------------------------------------------------------

    //check point location of points
    if (check_point_location(arr, plist))
    {
        std::cout << "ERROR in check_point_location."<<std::endl<<std::endl;
        return (false);
    }
    std::cout << std::endl;

    return (true);
}

/*!
  */
bool test4()
{

    //-Initialization-------------------------------------------------
    Arrangement_2 arr;
    Curve_list curve_list;
    Points_list plist;
    //-Curves definition----------------------------------------------
    /*
      xy=1,
      xy=-1,
      y=  (x+4)^2-2
    */
    curve_list.push_back(construct(P1,x));

    curve_list.push_back(construct(minusP1,x));

    curve_list.push_back(construct((x+4)*(x+4)-2,P1));

    insert_curves(arr, curve_list);

    //Single isolated vertex in target face.
    //BUG FIXED:Failed for specified point generator or segmentation fault.
    Traits_2::Construct_point_2 construct_point = traits.construct_point_2_object();
    //Isolated Vertex, query point and landmark point same coordinates.
    CGAL::insert_point(arr, construct_point(CORE::BigInt(-4),CORE::BigInt(-4)));
    //Isolated vertex on target face
    CGAL::insert_point(arr, construct_point(CORE::BigInt(4.1),CORE::BigInt(4.1)));

    //-Point generation------------------------------------------------------------
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(4), CORE::BigInt(4)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(4.1), CORE::BigInt(4.1)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(4), CORE::BigInt(-4)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(-4), CORE::BigInt(4)));
    plist.push_back(traits.construct_point_2_object()(CORE::BigInt(-4), CORE::BigInt(-4)));
    //-----------------------------------------------------------------------------

    //check point location of points
    if (check_point_location(arr, plist))
    {
        std::cout << "ERROR in check_point_location."<<std::endl<<std::endl;
        return (false);
    }
    std::cout << std::endl;

    return (true);
}

//-END TESTS----------------------------------------------------------------------------


int main(/*int argc, char * argv[]*/)
{
    int success = 0;

    CGAL::set_pretty_mode(std::cout); // for nice printouts.

    success+= test0();
    success+= test1();
    success+= test3();
    success+= test4();

    return success;

}

