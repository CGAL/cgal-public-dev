#include <CGAL/basic.h>

#ifdef CGAL_USE_GMP

  #include <CGAL/Gmpq.h>

  typedef CGAL::Gmpq                                    Number_type;

#else
  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>

  typedef CGAL::Quotient<CGAL::MP_Float>                Number_type;

#endif

#include <CGAL/Cartesian.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Timer.h>
#include <cstdlib>

#include <CGAL/Arr_segment_traits_2.h>
#include "Segment_reader.h"

#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_simple_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_point_location/Arr_landmarks_point_location_no_construction.h>
#include <CGAL/Arr_point_location/Arr_lm_random_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_grid_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_halton_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_middle_edges_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_specified_points_generator.h>
//#include <CGAL/Arr_triangulation_point_location.h>

typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits_2;

typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Curve_2                               Curve_2;
typedef Traits_2::X_monotone_curve_2                    Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;
typedef Arrangement_2::Halfedge_handle                  Halfedge_handle;
typedef Arrangement_2::Edge_const_iterator              Edge_const_iterator;
typedef Arrangement_2::Vertex_const_iterator            Vertex_const_iterator;

typedef CGAL::Arr_naive_point_location<Arrangement_2>     
                                                    Naive_point_location;
typedef CGAL::Arr_simple_point_location<Arrangement_2>     
                                                    Simple_point_location;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> 
                                                    Walk_point_location;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2> 
                                                    Lm_point_location;
typedef CGAL::Arr_random_landmarks_generator<Arrangement_2>
                                                    Random_lm_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Random_lm_generator> 
                                                    Lm_random_point_location;
typedef CGAL::Arr_grid_landmarks_generator<Arrangement_2>
                                                    Grid_lm_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Grid_lm_generator> 
                                                    Lm_grid_point_location;
typedef CGAL::Arr_halton_landmarks_generator<Arrangement_2>
                                                    Halton_lm_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Halton_lm_generator> 
                                                    Lm_halton_point_location;
typedef CGAL::Arr_middle_edges_landmarks_generator<Arrangement_2>
                                                    Middle_edges_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Middle_edges_generator> 
                                                    Lm_middle_edges_point_location;
typedef CGAL::Arr_landmarks_specified_points_generator<Arrangement_2>
                                                    Specified_points_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Specified_points_generator> 
                                                    Lm_specified_points_point_location;


typedef CGAL::Arr_landmarks_point_location_no_construction<Arrangement_2> 
                                                    Lm_point_location_no_c;
typedef CGAL::Arr_landmarks_point_location_no_construction<Arrangement_2, Random_lm_generator> 
                                                    Lm_random_point_location_no_c;
typedef CGAL::Arr_landmarks_point_location_no_construction<Arrangement_2, Grid_lm_generator> 
                                                    Lm_grid_point_location_no_c;
typedef CGAL::Arr_landmarks_point_location_no_construction<Arrangement_2, Halton_lm_generator> 
                                                    Lm_halton_point_location_no_c;
typedef CGAL::Arr_landmarks_point_location_no_construction<Arrangement_2, Middle_edges_generator> 
                                                    Lm_middle_edges_point_location_no_c;
typedef CGAL::Arr_landmarks_point_location_no_construction<Arrangement_2, Specified_points_generator> 
                                                    Lm_specified_points_point_location_no_c;




//typedef CGAL::Arr_triangulation_point_location<Arrangement_2> 
//                                                    Lm_triangulation_point_location;

// ===> Add new point location type here <===

typedef std::list<Point_2>                                Points_list;
typedef Points_list::iterator                             Point_iterator;
typedef std::list<Curve_2>                                Curve_list;
typedef std::vector<CGAL::Object>                         Objects_vector;
typedef Objects_vector::iterator                          Object_iterator;

// ===> Change the number of point-location startegies
//      when a new point location is added. <===
#define NUM_OF_POINT_LOCATION_STRATEGIES 15

template <class PL, class Out>
double run_pl (Points_list &plist, PL &pl, Out out) {
  CGAL::Timer timer;
  timer.reset(); 
  timer.start(); //START
  for (Point_iterator piter = plist.begin(); piter != plist.end(); piter++) {
    CGAL::Object obj = pl.locate (*piter);
    *out++ = obj;
  }
  timer.stop(); ///END
  return timer.time();
}

/*! */
int check_point_location (Arrangement_2 &arr, Points_list &plist)
{
  //init - all point locations
  CGAL::Timer            timer;

  Naive_point_location            naive_pl (arr);                 // 0

  Simple_point_location           simple_pl (arr);                // 1
  Walk_point_location             walk_pl (arr);                  // 2

  timer.reset(); timer.start();
  Lm_point_location               lm_pl (arr);                    // 3
  Lm_point_location_no_c          lm_noc_pl (arr);                // 4
  timer.stop(); 
  std::cout << "Lm (vert) construction took " << timer.time() <<std::endl;

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
  std::cout << "Grid lm construction took " << timer.time() <<std::endl;

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

  Specified_points_generator::Points_set points;
  //points.push_back(Point_2(1,1));
  //points.push_back(Point_2(2,2));
  timer.reset(); timer.start();
  //Specified_points_generator                specified_points_g(arr,points);
  Specified_points_generator                specified_points_g(arr);
  Lm_specified_points_point_location        specified_points_lm_pl (arr, &specified_points_g);  // 13
  Lm_specified_points_point_location_no_c   specified_points_lm_noc_pl (arr, &specified_points_g);  // 14
  timer.stop(); 
  std::cout << "Specified_points lm construction took " << timer.time() <<std::endl;

/*
  timer.reset(); timer.start();
  Lm_triangulation_point_location        triangulation_lm_pl (arr);  // 9
  timer.stop(); 
  std::cout << "Triangulation lm construction took " << timer.time() <<std::endl;
*/
  std::cout << std::endl;

  // ===> Add new point location instance here. <===

  Objects_vector                  objs[NUM_OF_POINT_LOCATION_STRATEGIES];
  Object_iterator                 ob_iter[NUM_OF_POINT_LOCATION_STRATEGIES];
  Arrangement_2::Vertex_const_handle    vh_ref, vh_curr;
  Arrangement_2::Halfedge_const_handle  hh_ref, hh_curr;
  Arrangement_2::Face_const_handle      fh_ref, fh_curr;
  
  Point_2               q;

  //LOCATE the points in the list using all PL strategies

  //std::cout << "Time in seconds" <<std::endl; ;
  std::cout << std::endl;
  
  std::cout << "Naive location took " << run_pl(plist, naive_pl, std::back_inserter(objs[0])) << std::endl;
  std::cout << "Simple location took " << run_pl(plist, simple_pl, std::back_inserter(objs[1])) << std::endl;
  std::cout << "Walk location took " << run_pl(plist, walk_pl, std::back_inserter(objs[2])) << std::endl;
  std::cout << "Landmarks (vertices) location took "  << run_pl(plist, lm_pl, std::back_inserter(objs[3])) << std::endl;
  std::cout << "Landmarks (vertices) location no construction took "  << run_pl(plist, lm_noc_pl, std::back_inserter(objs[4])) << std::endl;
  std::cout << "Random LM location took "  << run_pl(plist, random_lm_pl, std::back_inserter(objs[5])) << std::endl;
  std::cout << "Random LM location no construction took "  << run_pl(plist, random_lm_noc_pl, std::back_inserter(objs[6])) << std::endl;
  std::cout << "Grid LM location took "  << run_pl(plist, grid_lm_pl, std::back_inserter(objs[7])) << std::endl;
  std::cout << "Grid LM location no construction took "  << run_pl(plist, grid_lm_noc_pl, std::back_inserter(objs[8])) << std::endl;
  std::cout << "Halton LM location took "  << run_pl(plist, halton_lm_pl, std::back_inserter(objs[9])) << std::endl;
  std::cout << "Halton LM location no construction took "  << run_pl(plist, halton_lm_noc_pl, std::back_inserter(objs[10])) << std::endl;
  std::cout << "Middle edges LM location took "  << run_pl(plist, middle_edges_lm_pl, std::back_inserter(objs[11])) << std::endl;
  std::cout << "Middle edges LM location no construction took "  << run_pl(plist, middle_edges_lm_noc_pl, std::back_inserter(objs[12])) << std::endl;
  std::cout << "Specified points LM location took " << run_pl(plist, specified_points_lm_pl, std::back_inserter(objs[13])) << std::endl;
  std::cout << "Specified points LM location no construction took " << run_pl(plist, specified_points_lm_noc_pl, std::back_inserter(objs[14])) << std::endl;


/*
  timer.reset(); 
  timer.start(); //START
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = triangulation_lm_pl.locate (q);
    objs[8].push_back(obj);
  }
  timer.stop(); ///END
  std::cout << "Triangulation LM location took " << timer.time() <<std::endl;
*/
  std::cout << std::endl;

  // ===> Add a call to operate the the new point location. <===

  //END LOCATION 
  int pls_num = NUM_OF_POINT_LOCATION_STRATEGIES;
  int pl_index;
  int result = 0;

  //init all obejct iterators
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

  //assign and check results
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
            std::cout << ", an unknowen object returned instead of a face"<< std::endl;
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
      //if (fh_ref->is_unbounded())
      //  std::cout << "Unbounded face." << std::endl;
      //else
      //  std::cout << "Face." << std::endl;
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
            std::cout << ", an unknowen object returned instead of an halfedge"<< std::endl;
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

/*! */
int read_points(const char * points_filename, Points_list &plist)
{
  //read points from file into list
  std::ifstream inp_pnt_file(points_filename);
  if (!inp_pnt_file.is_open()) 
  {
    std::cerr << "Cannot open file " << points_filename << "!" << std::endl;
    return (-1);
  }

  int points_count = 0;
  inp_pnt_file >> points_count;

  for (int i = 0; i < points_count; i++) {
    Number_type x, y;
    inp_pnt_file >> x >> y;
    Point_2 pnt(x, y);
    plist.push_back(pnt);
  }

  return 0;
}

bool test(const char* curves_filename, const char* points_filename)
{
  //read curves and insert them into the arrangement 
  Segment_reader<Traits_2> reader;

  CGAL::Bbox_2             bbox;
  Curve_list               curve_list;
  reader.read_data(curves_filename, std::back_inserter(curve_list), bbox);  

  //insert all curves into the arrangement
  CGAL::Timer timer;
  timer.reset(); 
  timer.start(); //START
  Arrangement_2  arr;
  insert (arr, curve_list.begin(), curve_list.end());
  timer.stop(); ///END
  std::cout << "Arrangement aggregate construction took " 
            << timer.time() <<std::endl;  
  // Print the size of the arrangement.
  std::cout << "V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;

  //read point and insert them into a list of points
  Points_list           plist;

  //read points from file into list
  if (read_points(points_filename, plist))
  {
    std::cout << "ERROR in read_points."<<std::endl<<std::endl;
    return (false);
  }

  //check point location of points
  if (check_point_location(arr, plist))
  {
    std::cout << "ERROR in check_point_location."<<std::endl<<std::endl;
    return (false);
  }
  std::cout << std::endl;
  return (true);
}

int main (int argc, char * argv[])
{
   //get arguments
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] <<" curve_file pnt_file" << std::endl;
    std::cout << "curve_file  - the input curves file" << std::endl;
    std::cout << "pnt_file    - the input query points" << std::endl;
    std::exit(-1);
  }
  int success = 0;
  for(int i=1; i<argc; i+=2)
  {
    const char * curves_filename = argv[i];
    const char * points_filename = argv[i+1];

    if(!test(curves_filename, points_filename))
    {
      std::cout<<"ERROR : "<<argv[0]<<" "<<argv[i]<<" "<<argv[i+1]<<std::endl;
      success = -1;
    }
  }

  return (success);
}
