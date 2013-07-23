#include <CGAL/basic.h> 

#include <SV/Swept_volume_with_vhull_3.h>
#include <SV/Mesh_domain_3.h>

#include <Swept_volume_3_mesh_function.h>
#include <Mesh_parameters.h> 
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <iostream> 
#include <fstream> 

#include <CGAL/Real_timer.h>
#include <SV/io.h>

#include <time.h>

#if USE_OMP
#include <omp.h>
#endif 

//#undef CGAL_SURFACE_MESHER_PROFILE

// test swept volume meshing scheme as used in the demo. 
int main(int argc, char **argv){
  {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel    Kernel; 
    typedef SV::Swept_volume_with_vhull_3<Kernel>                  Volume; 
    typedef SV::Mesh_domain_3<Volume>                              Mesh_domain; 
    typedef Swept_volume_3_mesh_function<Mesh_domain>              Mesh_function; 
    typedef Mesh_function::C3t3                                    C3t3; 
  
	  typedef typename Kernel::Point_3 Point_3;
    typedef typename Kernel::Aff_transformation_3 AT_3;

    CGAL::Real_timer total; 
    total.start(); 

///


#if 0
  typedef Kernel::Point_3 Point_3;
	//typedef CGAL::Aff_transformation_3<Kernel> AT3;
	typedef Kernel::Aff_transformation_3 AT_3;

	Point_3 a(1,0,3);
	Point_3 b(1,2,3);

	AT_3 A(1,0,0,1,
				0,1,0,1,
				0,0,1,1,
				      1);
	a = A.transform(a);
	std::cout << a.x() << " " << a.y() << " " << a.z() << std::endl;	

		




#else

#if 1
    std::string tarr_file = "../../data/Swept_volume/bunny/track.tarr";
    std::string off_file  = "../../data/Swept_volume/bunny/bunny.off";
#endif 

    // cmake -DCMAKE_BUILD_TYPE=Debug

#if 1
    int resolution_D = 9;
    if (argv[1] != NULL)
      off_file = argv[1];
    else
    {
      printf("Usage: SVTool <path to generator (.off)> <path to track (.tarr)> <resolution level (9,...,13)> <downstep level>\n");
      printf("Example: SVTool ../../engine.off ../../Tracks/vibration.path 12 \n");
      printf("OR use downstep level if memory becomes scarce\n");
      printf("Example: SVTool ../../engine.off ../../Tracks/vibration.path 12 1\n");
      exit(0);

    }
    if (argv[2] != NULL)
      tarr_file = argv[2];
    if (argv[3] != NULL)
      resolution_D = atoi(argv[3]);
    int downstep = 0; 
    if (argv[4] != NULL)
      downstep= atoi(argv[4]);
#endif

#if USE_OMP
    int threads = omp_get_max_threads();
#else
    int threads = 1; 
#endif 

    printf("Number of threads: %d \n",threads);

#ifdef NDEBUG    
    Volume volume(tarr_file, off_file, resolution_D, downstep, threads); 
#else
    Volume volume(tarr_file, off_file, resolution_D, downstep, threads);
#endif 
    std::cerr << "Dateien geladen!"<< std::endl; 

#if 1
    Mesh_parameters param;
    param.facet_angle = 0;
    param.facet_sizing = 0;
    param.facet_approx = 0;
    param.tet_sizing = 0;
    param.tet_shape = 0;

    C3t3 c3t3; 


    // will be deleted by the destructor of Volume_mesh_function which is not nice. ITS A BUG ! 
    Mesh_domain* domain = new Mesh_domain(volume);
    Mesh_function mesh_function(c3t3,domain,param);

    // could be stopped by a second thread, e.g. depending on a timer using mesh_function.stop();
    mesh_function.launch(-1);

    Point_3 e1 = volume.getBackTrafo().transform(Point_3(1,0,0));
    Point_3 e2 = volume.getBackTrafo().transform(Point_3(0,1,0));
    Point_3 e3 = volume.getBackTrafo().transform(Point_3(0,0,1));
    //std::cerr << e1 << " " << e2 << " " << e3 << std::endl;

    std::ofstream os("result_for_swept_volume_with_vhull_3.off");
    SV::save_as_off(c3t3,os,volume.getBackTrafo()); 

    total.stop(); 
    std::cout << "Result exported to: ./result_for_swept_volume_with_vhull_3.off" << std::endl;
    std::cout << "Mesh number of facets : " << c3t3.number_of_facets() << std::endl;
    std::cout << "Total time: " << total.time() << std::endl;


    // delete domain; 
#endif 
  }
  std::cerr << " Final memory usage     " << SV::memory_usage()*100 << " %" << std::endl ;
  std::cerr << " Final memory usage max " << SV::memory_usage_max()*100 << " %" << std::endl ;

  return 0; 

#endif
}
