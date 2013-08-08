#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Swept_volume_domain_3.h>

#include <iostream> 
#include <fstream> 


// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Swept_volume_domain_3< K > Mesh_domain; 
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;


int main(int argc, char **argv){
  {
    typedef typename K::Point_3 Point_3;
    typedef typename K::Aff_transformation_3 AT_3;

    std::string trj_file  = "track.trj";
    std::string off_file  = "bunny.off";
    
    std::vector<Point_3> vertices;
    std::vector< CGAL::cpp11::tuple<int, int, int> > indices; 
    std::vector<Aff_transformation_3> trajectory;
    
    // load Off File
    std::ifstream inStream(off_file);
    std::string fileHeader;
    std::getline(inStream, fileHeader);

    int N_triangles, N_vertices, dummy; 
    inStream >> N_vertices >> N_triangles >> dummy; 

    for(int i = 0; i < N_vertices; i++){
      double x, y, z;
      inStream >> x >> y >> z;
      vertices.push_back(Point_3(x, y, z));
    }

    for (int i = 0; i < N_triangles; i++){
        int x, y, z;
        inStream >> dummy >> x >> y >> z; 
        indices.push_back(std::make_tuple(x,y,z));
    }

    // load trj file
    std::ifstream inStream(trj_file);

    int N_trajec;
    inStream >> N_trajec;
    for(int i = 0; i < N_trajec; i++) {
      double t0,t1,t2,m00,m01,m02,m10,m11,m12,m20,m21,m22;
      inStream >> m00 >> m10 >> m20
               >> m01 >> m11 >> m21
               >> m02 >> m12 >> m22
               >>  t0 >>  t1 >>  t2;
      
      trajectory.push_back( AT_3( 
                                m00, m01, m02, t0,
                                m10, m11, m12, t1,
                                m20, m21, m22, t2,
                                1.0 ));
    }
    
    Swept_volume_domain_3 sv_domain(vertices, indices, trajectory, 0.005, false);
    
    
    CGAL::Mesh_criteria_3<Tr> normal_mesh_criteria( 
        facet_angle=30, facet_size=5, facet_distance=1.5,
        cell_radius_edge_ratio=0, cell_size=0);
    
    Swept_volume_domain_3::Swept_volume_criteria_3< MeshCriteria_3 >
      mesh_criteria_with_guarantees = 
      sv_domain.swept_volume_criteria_3_object(normal_mesh_criteria); 
    
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
      
    std::ofstream os("result_for_swept_volume_with_vhull_3.off");
    c3t3.output_boundary_to_off(os);
  }
  
