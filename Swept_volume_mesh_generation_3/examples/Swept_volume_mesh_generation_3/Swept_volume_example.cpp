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


    std::string trj_file  = "track.trj";
    std::string off_file  = "bunny.off";
    int resolution_D = 9;
    double epsilon = 0.005;
    int downstep = 0; 
    std::vector<Point_3> vertices;
    std::vector< CGAL::cpp11::tuple<int, int, int> > indices; 
    std::vector<Aff_transformation_3> trajectory;
    
    double epsilon; 
    bool downstep = false;

    // load Off File
    int N_triangles = 0;
    int N_vertices = 0;

    std::ifstream inStream(off_file);
    std::string fileHeader;
    std::getline(inStream, fileHeader);

    inStream >> N_vertices;
    inStream >> N_triangles;
    int tmp;
    inStream >> tmp;

    for(int i = 0; i < N_vertices; i++) 
    {
        float x, y, z;
        inStream >> x;
        inStream >> y;
        inStream >> z;
        vertices.push_back(Point_3(x, y, z));
    }

  for (int i = 0; i < N_triangles; i++)
    {
        int x, y, z;
        inStream >> tmp; //3
        inStream >> x;
        inStream >> y;
        inStream >> z;
        indices.push_back(std::make_tuple(x,y,z));
    }

    // load trj file

    std::ifstream inStream(trj_file);

    inStream >> N_trajec;
    for(int i = 0; i < N_trajec; i++) 
    {
        double t0=0,t1=0,t2=0;
        double m00=0,m01=0,m02=0,m10=0,m11=0,m12=0,m20=0,m21=0,m22=0;
        inStream >> m00; inStream >> m10; inStream >> m20;
        inStream >> m01; inStream >> m11; inStream >> m21;
        inStream >> m02; inStream >> m12; inStream >> m22;
        inStream >>  t0; inStream >>  t1; inStream >>  t2;

        trajectory.push_back(AT_3( 
            m00, m01, m02, t0,
            m10, m11, m12, t1,
            m20, m21, m22, t2,
                          1.0 ))  ;
    }
   

// TODO Swept_volume_criteria
    Mesh_parameters param;
    param.facet_angle = 0;
    param.facet_sizing = 0;
    param.facet_approx = 0;
    param.tet_sizing = 0;
    param.tet_shape = 0;

    C3t3 c3t3; 

    Swept_volume_domain_3 sv_domain(vertices, indices, trajectory, 0.005, false);

    Mesh_function mesh_function(c3t3,_sv_domain,param);

    // could be stopped by a second thread, e.g. depending on a timer using mesh_function.stop();
    mesh_function.launch(-1);

    std::ofstream os("result_for_swept_volume_with_vhull_3.off");
    SV::save_as_off(c3t3,os); 

    total.stop(); 

}
