#include <CGAL/Cartesian_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Timer.h>
#include <time.h> 
#include <fstream>

typedef CGAL::Cartesian_d<double>                         K;
typedef K::Point_d                                        Point_d;
typedef CGAL::Search_traits_d<K>                          Traits;
typedef CGAL::Random_points_in_cube_d<Point_d>            Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator>   N_Random_points_iterator;
typedef CGAL::Kd_tree<Traits>                             Tree;
typedef CGAL::Fuzzy_iso_box<Traits>                       Fuzzy_iso_box;
typedef std::list<Point_d>                                Points_list;

int main(int argc, char* argv[]) {
  const int D = 3;
  const int N = 10000000;

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

  Points_list result;

  Timer_overall.start();
  
  for(int bound = 100; bound<2000; bound+=100)
  {
    Timer_iteration.start();

    //creates a random point iterator. Random numbers ranging between 0 - bound*10 
    //bound also shows the bounding box of the Kd-tree
    Random_points_iterator rpit(3, (float)(bound*10) );

    // Insert N points in the tree
    Timer_tree_build_iteration.start();
    Timer_tree_build_overall.start();
    
    //this constructor call of the Kd_tree only initializes the points and does not build the Kd_tree.
    Tree tree(N_Random_points_iterator(rpit,0), N_Random_points_iterator(rpit, (bound*bound) ) );
    //Builds the kd tree with the random points initialized.
    tree.build();
    
    Timer_tree_build_iteration.stop();
    Timer_tree_build_overall.stop();

    //searching Begins
    Timer_search_iteration.start();
    Timer_search_overall.start();

    for(int bounding_box_width=10; bounding_box_width<bound; bounding_box_width+=2)
    {
      // define range query objects
      double  pcoord[D] = { 0.0, 0.0, 0.0 };
      double  qcoord[D] = { (float) (bounding_box_width), (float) (bounding_box_width), (float) (bounding_box_width)};
      Point_d p(D, pcoord, pcoord+D);
      Point_d q(D, qcoord, qcoord+D);
      
      //create a fuzzy box
      Fuzzy_iso_box fib(p, q);

      //Search function of Kd_tree will build the tree if not only built and the perform the search. 
      tree.search(std::back_inserter(result), fib);

      if(Verbose)
      {
        // //Print the putput to Standard Output.
        std::cout << "points approximately in fuzzy range query ";
        std::cout << "[<0.0,0.0, 0.0>,<" << bounding_box_width << "," << bounding_box_width << "," << bounding_box_width<< ">] are:" << std::endl;
        std::copy (result.begin(),result.end(), std::ostream_iterator<Point_d>(std::cout,"\n") );
        std::cout << std::endl;
      }
    }
    Timer_search_iteration.stop();
    Timer_search_overall.stop();

    result.clear();

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
    
    //reset timers
    Timer_iteration.reset();
    Timer_tree_build_iteration.reset();
    Timer_search_iteration.reset();
  } //forloop bound

  Timer_overall.stop();
  
  std::cerr << "\n\nOverall Tree Building took: " << Timer_tree_build_overall.time() << " sec" <<std::endl;
  std::cerr << "Overall Tree Searching took: " << Timer_search_overall.time() << " sec" <<std::endl;
  std::cerr << "Overall the process took " << Timer_overall.time() << " sec." << std::endl;

  if(Write_to_file)
  {
    Outfile << "\n\nOverall Tree Building took: " << Timer_tree_build_overall.time() << " sec" <<std::endl;
    Outfile << "Overall Tree Searching took: " << Timer_search_overall.time() << " sec" <<std::endl;
    Outfile << "Overall the process took " << Timer_overall.time() << " sec." << std::endl;
    Outfile.close();
  }
  return 0;
}