#include <CGAL/Cartesian.h>

#include <iostream>
#include <iterator>
#include <ctime>
#include <cassert>
#include <cstdlib>
#include <list>

#include <CGAL/kdtree_d.h>
#include <CGAL/Timer.h>
#include <time.h> 
#include <fstream>

template <int  DIM>
class Point_float_d
{
private:
  double   vec[ DIM ];

public:
  Point_float_d()
  {
    for  ( int ind = 0; ind < DIM; ind++ )
      vec[ ind ] = 0;
  }

  int dimension() const
  {
    return  DIM;
  }

//not essential by specification but needed for initializing a general d-point
  void set_coord(int k, double x)
  {
    assert( 0 <= k  &&  k < DIM );
    vec[ k ] = x;
  }

  double  & operator[](int k)
  {
    assert( 0 <= k  &&  k < DIM );
    return  vec[ k ];
  }

  double  operator[](int k) const
  {
    assert( 0 <= k  &&  k < DIM );
    return  vec[ k ];
  }
};

// not essential by specification but nice to have
template <int DIM>
std::ostream &operator<<(std::ostream &os, const Point_float_d<DIM> &p)
{
  std::cout << "(";
  for(int i = 0; i < DIM; i++)
    {
      std::cout << p[i] ;
      if (i < p.dimension() - 1) std::cout << ", ";
    }
  std::cout << ")";
  return os;
}

typedef Point_float_d<4>  point;
typedef CGAL::Kdtree_interface<point>   kd_interface;
typedef CGAL::Kdtree_d<kd_interface>    kd_tree;
typedef kd_tree::Box  box;
typedef std::list<point>                points_list;

//RANDOM FUNCTIONS
// dblRand - a random number between 0..1
#ifndef  RAND_MAX
#define  RAND_MAX    0x7fffffff
#endif

inline double dblRand( void )
{
    return  (double)std::rand() / (double)RAND_MAX;
}

void random_points( int  num, points_list &l, int DIM, double bound)
{
  double  x;

  for  (int j = 0;  j < num; j++)
    {
      point p;
      for (int i=0; i<=DIM; i++)
        {
          x = dblRand()*bound;
          p.set_coord(i,x);
        }
      l.push_front(p);
    }
}

int main(int argc, char* argv[])
{

   //flags that will be enabled by the command line arguments.
  bool Write_to_file = false, Verbose = false;

  std::ofstream Outfile;

  if( argc<=3 )
  {
    if( argc <= 3 && argc!=1) //file name specified
    {
      Write_to_file = true;
      std::string File_name = argv[1];

      Outfile.open(File_name);
      if(!Outfile.is_open())
      {
        std::cerr << "Can not open the file " << argv[1] << std::endl;
        return -1;
      }
      else
      {
        time_t rawtime;
        struct tm * timeinfo;

        time (&rawtime);
        timeinfo = localtime (&rawtime);

        Outfile << "Benchmark Started at " << asctime(timeinfo) << std::endl;
      }
    }

    if(argc == 3) // activate verbose
    {
      std::string Arg = argv[2];
      if( Arg.compare("-v") ==0 )
        Verbose = true;
      else
      {
        std::cerr << "Unknown Command line arguments. First argument=Output file name(optional) and Second_argument=-v (optinal)" << std::endl;
        return -1;
      }
    }
  }
  else
  {
    std::cerr << "Unknown Command line arguments. First argument=Output file name(optional) and Second_argument=-v (optinal)" << std::endl;
  }

  CGAL::Timer Timer_iteration, Timer_overall,
              Timer_tree_build_iteration, Timer_tree_build_overall, 
              Timer_search_iteration, Timer_search_overall;


  Timer_overall.start();

  CGAL::Kdtree_d<kd_interface>  tree(3);

  std::srand( (unsigned)time(NULL) );

   points_list  l , res;

  for(int bound = 100; bound<2000; bound+=100)
  {
    // N=bound*bound  container='l' Dimension=3   Random number limit=bound*10
    random_points( bound*bound, l , 3, bound*10 );

    Timer_iteration.start();
    
    // Building the tree for the random points
    Timer_tree_build_iteration.start();
    Timer_tree_build_overall.start();
    
    tree.build( l );
    
    Timer_tree_build_iteration.stop();
    Timer_tree_build_overall.stop();

    // Checking validity
    if  ( ! tree.is_valid() )
      tree.dump();
    assert( tree.is_valid() );

    // Searching the box r
    point p,q;
    
    Timer_search_iteration.start();
    Timer_search_overall.start();
    for(int bounding_box_width=10; bounding_box_width<bound; bounding_box_width+=2)
    {
      for (int k=0; k<=3; k++)
      {
        p.set_coord(k,0);
        q.set_coord(k,bounding_box_width);
      }
    
      box r(p, q, 3);
      tree.search( std::back_inserter(res), r );

      if(Verbose)
      {
        std::cout << "[<0.0,0.0, 0.0>,<" << bounding_box_width << "," << bounding_box_width << "," << bounding_box_width<< ">] are:" << std::endl;
        std::copy (res.begin(),res.end(), std::ostream_iterator<point>(std::cout,"\n") );
        std::cout << std::endl;
      }
      res.clear();
    }
    Timer_search_iteration.stop();
    Timer_search_overall.stop();

    tree.delete_all();

    Timer_iteration.stop();

    std::cerr << "Tree with Limit: "<< bound*10 <<" took " <<Timer_iteration.time() <<" sec." 
              <<"\t\tTree Build took: " <<Timer_tree_build_iteration.time() << " sec." 
              <<"\t\tTree Search took: " <<Timer_search_iteration.time() << " sec."<<std::endl;

    if(Write_to_file)
    {
      Outfile << "Tree with Limit: "<< bound*10 <<" took " <<Timer_iteration.time() <<" sec." 
              <<"\t\tTree Build took: " <<Timer_tree_build_iteration.time() << " sec." 
              <<"\t\tTree Search took: " <<Timer_search_iteration.time() << " sec."<<std::endl;
    }
    
    Timer_iteration.reset();
    Timer_tree_build_iteration.reset();
    Timer_search_iteration.reset();
    
    //clear the lists
    l.clear();
    res.clear();
  } 

  Timer_overall.stop();

  std::cerr << "\n\nOverall Tree Building took: " << Timer_tree_build_overall.time() << " sec" <<std::endl;
  std::cerr << "Overall Tree Searching took: " << Timer_search_overall.time() << " sec" <<std::endl;
  std::cerr << "Overall the process took " << Timer_overall.time() << " sec." << std::endl;

  if(Write_to_file)
  {
    Outfile << "\n\nOverall Tree Building took: " << Timer_tree_build_overall.time() << " sec" <<std::endl;
    Outfile << "Overall Tree Searching took: " << Timer_search_overall.time() << " sec" <<std::endl;
    Outfile << "Overall the process took " << Timer_overall.time() << " sec." << std::endl;
  }

  return 0;
}
