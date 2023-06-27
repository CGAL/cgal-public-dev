#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Variational_shape_reconstruction.h>
#include <iostream>

#include <CGAL/Point_set_3/IO/XYZ.h>
// CGAL
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Point_set_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;


typedef Kernel::FT                  FT;
typedef Kernel::Point_3             Point;
typedef Kernel::Vector_3            Vector;
typedef Kernel::Triangle_3          Triangle;
typedef Kernel::Plane_3             Plane;

typedef Kernel::Point_2             Point_2;
typedef Kernel::Segment_2           Segment_2;

typedef CGAL::First_of_pair_property_map<std::pair<Point, Vector>>                     Point_map;
typedef CGAL::Second_of_pair_property_map<std::pair<Point, Vector>>                    Normal_map;
typedef CGAL::Point_set_3< Point, Vector > Pointset;
Pointset load(std::string m_filename)
{

        std::ifstream fstream(m_filename);


		Pointset pointset;

        if (!fstream || !CGAL::IO::read_XYZ( fstream,pointset))
        {
            std::cerr << "Error: cannot read file " << m_filename << std::endl;
        } 

		return  pointset;
}
int main(int argc, char** argv)
{	

    const std::string fname = argc > 1 ? argv[1] :  CGAL::data_file_path("../data/piece_meca.xyz");
    Pointset pointset;
    if (!CGAL::IO::read_XYZ( fname,pointset))
    {
        std::cerr << "Error: cannot read file " << std::endl;
        return 0;
    } 
    const size_t generators = 30;
    const size_t steps = 10;
    const double split_threshold =0.01;

    // reconstruction
    const double  dist_ratio = 0.001;
	const double  fitting = 0.43;
	const double  coverage = 0.27;
	const double  complexity = 0.3;

    size_t new_generators = generators; 
    size_t iteration = 0 ;

    
	 for( Pointset::const_iterator pointset_it = pointset.begin(); pointset_it != pointset.end(); ++ pointset_it )
        {
            const auto point = pointset.normal(*pointset_it);
        }
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	qem::Variational_shape_reconstruction manager(pointset,generators);
        std::chrono::steady_clock::time_point begin_clustering = std::chrono::steady_clock::now();
		//manager.region_growing(steps);
    while(new_generators > 5 )
    {
        
		manager.region_growing(steps);
        new_generators = manager.guided_split_clusters(split_threshold, iteration++);
		Pointset point_cloud = manager.get_point_cloud_clustered();
		std::ofstream edge_file;
		edge_file.open("clustering_"+std::to_string(iteration)+".ply");

		edge_file << "ply\n"
					<< "format ascii 1.0\n"
					<< "element vertex " << point_cloud.size() << "\n"
					<< "property float x\n"
					<< "property float y\n"
					<< "property float z\n"
					<< "property uchar red\n"
					<< "property uchar green\n"
					<< "property uchar blue\n"
					<< "end_header\n";

		for(Pointset::const_iterator it = point_cloud.begin(); it != point_cloud.end(); ++ it)
		{
			auto point = point_cloud.point(*it);
			edge_file << point.x() << " " << point.y() << " " << point.z() << " ";
			auto normal = point_cloud.normal(*it);
			edge_file << static_cast<int>(255*normal.x()) << " " << static_cast<int>(255*normal.y()) << " " << static_cast<int>(255*normal.z()) << std::endl;
		}
		edge_file.close();
    }

    
    

    //CGAL::IO::write_XYZ("clustering_"+fname+".txt", point_cloud);



    
     std::chrono::steady_clock::time_point end_clustering = std::chrono::steady_clock::now();
    std::cerr << "Clustering " << std::chrono::duration_cast<std::chrono::microseconds>(end_clustering - begin_clustering).count()/1000 << "[ms]" << std::endl;
    
    // Reconstruction
    manager.reconstruction(dist_ratio, fitting, coverage, complexity);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cerr << "Algo " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000 << "[ms]" << std::endl;
	Polyhedron mesh = manager.get_reconstructed();

        std::ofstream mesh_file;
        mesh_file.open("mesh.off");
        CGAL::write_off(mesh_file, mesh);
        mesh_file.close();
    /*
   
    auto points = load(fname);
	//save(points);
	qem::Variational_shape_reconstruction manager2;
	manager2.read_xyz(fname.c_str());
	manager2.initialize_poles(generators);
	//manager2.read_xyz(""../data/piece_meca.xyz");
	//manager2.init_random_poles(30);
	for(int p = 0 ; p < 6;p++)
	{
		for(int i = 0 ; i < 10;i++)
		{
			manager2.region_growing(true);
			auto flag = manager2.update_poles();

			if(!flag)
				break;
		}
		//Clustering();
		manager2.guided_split_clusters(0.01);
	}
	// adjacent
	std::vector<float> adjacent_edges;
	std::vector<float> candidate_facets;
	std::vector<float> candidate_normals;
	// A
	manager2.create_adjacent_edges();
	manager2.update_adjacent_edges(adjacent_edges);
	//B
	manager2.create_candidate_facets();
	manager2.update_candidate_facets( candidate_facets, candidate_normals);
	//C

	manager2.mlp_reconstruction(dist_ratio, fitting, coverage, complexity);
	std::vector<float> fit_facets;
	std::vector<float> fit_normals;
	std::vector<float> fit_soup_facets;
	std::vector<float> fit_soup_normals;
	manager2.update_fit_surface(fit_facets, fit_normals);
	manager2.update_fit_soup(fit_soup_facets, fit_soup_normals);
	*/
	
    
	return 0;
}


