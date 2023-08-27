#include "io.h"
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Splitters.h>

#include "trianglefit.h"

#include "clustering.h"

namespace qem {
typedef typename std::unordered_set<std::pair<int, int>, HashPairIndex, EqualPairIndex>     IntPairSet; 
typedef std::vector<IntList>                                    IntListList;
typedef CGAL::Bbox_3   Bbox;


class Variational_shape_reconstruction
{
    private:

        //generators
        size_t m_generator_count;
        std::vector<int> m_generators;

        // geometry
        std::vector<std::pair<Point, size_t>> m_points;
        Pointset pointset_;

        // knntree
        KNNTree m_tree;
        int m_num_knn = 12;
        

        //Region growing
        std::map<int, int> m_vlabels;
        std::vector<QEM_metric> m_generators_qem;

        //init
        Bbox            m_bbox;
        double          m_diag;
        double          m_spacing;

        TriangleFit              m_triangle_fit;
        std::vector<int> m_generators_count;
        std::shared_ptr<Clustering> m_cluster;
        int m_verbose_level;


    public:
    Variational_shape_reconstruction(const Pointset& pointset,
    int generator_count,
    double euclidean_distance_weight,
    int verbose_level,
    int initialization_algorithm
    ) :
    m_verbose_level(verbose_level),
    m_generator_count(generator_count)
    {
        pointset_ = pointset;
        load_points(pointset_);
        compute_bounding_box();   

        // init kdtree
        m_tree.clear();
        m_tree.insert(m_points.begin(), m_points.end());     
                                                            
        // compute average spacing
        m_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(m_points, 6,
         CGAL::parameters::point_map(CGAL::First_of_pair_property_map<std::pair<Point, std::size_t>>()));            

        m_cluster = std::make_shared<Clustering>(pointset, m_num_knn,euclidean_distance_weight, verbose_level);

        std::cout<<"Initialization.\n";
        switch(initialization_algorithm)
        {
            case 1:
                init_random_generators_kmeanspp();
            break;
            case 2:
                init_random_generators_kmeans_farthest();
            break;
            default:
                init_random_generators();
            break;

        }
        m_cluster->initialize_qem_map(m_tree);
        m_cluster->initialize_vertex_qem(m_tree, m_generators);    
        if(m_verbose_level==3)
        {
            auto point_cloud = get_point_cloud_clustered();
            std::ofstream clustering_file;
            clustering_file.open("clustering_init.ply");
            clustering_file << "ply\n"
                        << "format ascii 1.0\n"
                        << "element vertex " << m_generators.size() << "\n"
                        << "property float x\n"
                        << "property float y\n"
                        << "property float z\n"
                        << "property uchar red\n"
                        << "property uchar green\n"
                        << "property uchar blue\n"
                        << "end_header\n";
            std::vector<Vector> colors;
            for(int i = 0 ; i < m_generators.size(); i++)
            {
                double r = (double) rand() / (RAND_MAX);
                double g = (double) rand() / (RAND_MAX);
                double b = (double) rand() / (RAND_MAX);
                colors.push_back(Vector(r,g,b));
            }
            int point_index =0;
            int generator_index=0;
            for(Pointset::const_iterator it = point_cloud.begin(); it != point_cloud.end(); ++ it)
            {
                if(std::find(m_generators.begin(),m_generators.end(),point_index)!=m_generators.end())
                {
                    auto point = point_cloud.point(*it);
                    clustering_file << point.x() << " " << point.y() << " " << point.z() << " ";
                    auto normal = colors[generator_index++];
                    int r = static_cast<int>(255*normal.x());
                    int g = static_cast<int>(255*normal.y());
                    int b = static_cast<int>(255*normal.z());
                    clustering_file << r << " " << g << " " << b << std::endl;
                }
                point_index++;
            }
            clustering_file.close();
            std::cout<<"clustering at initialization written to disk.\n";
        }

    }
    /// @brief load the points from a pointset to the m_points list
    /// @param pointset 
    void load_points(const Pointset& pointset)
    {
        for( Pointset::const_iterator pointset_it = pointset.begin(); pointset_it != pointset.end(); ++ pointset_it )
        {
            const auto point = pointset.point(*pointset_it);
            m_points.push_back(std::make_pair(point, std::distance<Pointset::const_iterator>(pointset.begin(), pointset_it)));    
        }
        if(m_verbose_level > 1)
            std::cout << "Number of points: " << pointset_.size() << std::endl;
    }
    /// @brief Compute the bounding box of the pointset
    void compute_bounding_box()
    {
        // find bounding box
        // todoquestion : boost bbox over poinset 
        boost::function<Point(std::pair<Point, size_t>&)> pwi_it_to_point_it = boost::bind(&std::pair<Point, size_t>::first, _1);
        m_bbox = CGAL::bbox_3(  boost::make_transform_iterator(m_points.begin(), pwi_it_to_point_it), 
                                boost::make_transform_iterator(m_points.end(), pwi_it_to_point_it));
        m_diag = std::sqrt(CGAL::squared_distance(Point(m_bbox.min(0), m_bbox.min(1), m_bbox.min(2)),
                                                Point(m_bbox.max(0), m_bbox.max(1), m_bbox.max(2))));
        if(m_verbose_level > 1)
            std::cout << "Diagonal of bounding box: " << m_diag << std::endl;
    }
    /// @brief Initialize with random generators
    void init_random_generators()
    {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        std::vector<int> num_range(pointset_.size());
        std::iota(num_range.begin(), num_range.end(), 0);

        std::random_device rd;
        std::mt19937 g(27);
        std::shuffle(num_range.begin(), num_range.end(), g);

        std::set<int> selected_indices;
        for(int i = 0; i < m_generator_count; i++)
            selected_indices.insert(num_range[i]);
            
        for(auto &elem: selected_indices) {
            m_generators.push_back(elem);
        }   
        if(m_verbose_level > 1)
            std::cout << "Number of random poles: " << m_generators.size() << std::endl;

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        if(m_verbose_level > 1)
            std::cerr << "Random poles in " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[us]" << std::endl;
        
    }
    /// @brief Initiatlize starting with one random generator and kmeans++
    void init_random_generators_kmeanspp()
    {
        std::vector<int> num_range(pointset_.size());
        std::iota(num_range.begin(), num_range.end(), 0);
        std::set<int> selected_indices;
        // one generator for k means ++
        for(int i = 0; i < 1; i++)
            selected_indices.insert(num_range[i]);
            
        for(auto &elem: selected_indices) {
            m_generators.push_back(elem);
        }   
        std::mt19937 gen;
        for(int i = 1; i < m_generator_count; i++)
        {
            std::vector<float> distance_list;
            for( Pointset::const_iterator pointset_it = pointset_.begin(); pointset_it != pointset_.end(); ++ pointset_it )
            {
                const auto point = pointset_.point(*pointset_it);
                float distance = std::numeric_limits<float>::max();
                // we iterate over previously selected generators
                for(int j = 0 ; j < m_generators.size(); j++)
                {
                    float distance_to_generator = CGAL::squared_distance(pointset_.point(m_generators[j]), point);
                    distance = std::min(distance, distance_to_generator);
                }
                distance_list.push_back(distance);
            }
            std::discrete_distribution<int> dist(std::begin(distance_list), std::end(distance_list));
            gen.seed(time(0));
            size_t index_max_distance = dist(gen);

            m_generators.push_back(index_max_distance);
        }

    }
    /// @brief Initiatlize starting with one random generator and kmeans farthest
    void init_random_generators_kmeans_farthest()
    {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        std::vector<int> num_range(pointset_.size());
        std::iota(num_range.begin(), num_range.end(), 0);

        std::random_device rd;
        std::mt19937 g(27);
        std::shuffle(num_range.begin(), num_range.end(), g);

        std::set<int> selected_indices;

        // one generator for k means farthest
        for(int i = 0; i < 1; i++)
            selected_indices.insert(num_range[i]);
            
        for(auto &elem: selected_indices) {
            m_generators.push_back(elem);
        }   

        for(int i = 1; i < m_generator_count; i++)
        {
            std::vector<float> distance_list;
            for( Pointset::const_iterator pointset_it = pointset_.begin(); pointset_it != pointset_.end(); ++ pointset_it )
            {
                const auto point = pointset_.point(*pointset_it);
                float distance = std::numeric_limits<float>::max();
                // we iterate over previously selected generators
                for(int j = 0 ; j < m_generators.size(); j++)
                {
                    float distance_to_generator = CGAL::squared_distance(pointset_.point(m_generators[j]),point);
                    distance = std::min(distance, distance_to_generator);
                }
                distance_list.push_back(distance);
            }
            auto max_distance_iterator = std::max_element(distance_list.begin(),distance_list.end());
            size_t index_max_distance = max_distance_iterator-distance_list.begin();
            m_generators.push_back(index_max_distance);
        }

        if(m_verbose_level > 1)
            std::cout << "Number of random poles: " << m_generators.size() << std::endl;

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        if(m_verbose_level > 1)
            std::cerr << "Random poles in " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[us]" << std::endl;
        
    }
    /// @brief region growing for n steps user defined
    /// @param steps 
    void region_growing(size_t steps)
    {
        bool flag = true;
        for(int i = 0 ; i < steps && flag; i++)
        {
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            m_vlabels.clear();
            m_generators_qem.clear();
            m_cluster->region_growing(m_vlabels, m_generators_qem, m_generators, true);
            assert(m_vlabels.size() == pointset_.size());

            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000;
            if(m_verbose_level > 1)
                std::cerr << "\nRegion growing in " << elapsed << "[ms]" << std::endl;
            flag = m_cluster->update_poles(m_vlabels,m_generators_qem, m_generators);
        }
    }
    /// @brief Split the cluster with the split_ratio
    /// @param split_ratio 
    /// @param iteration 
    /// @return total of added generators
    size_t guided_split_clusters(double split_ratio, size_t iteration) // batch splitting
    {
        std::cout << "Begin guided split..." << std::endl;
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        auto clusters_count = m_cluster->guided_split_clusters(m_vlabels,m_generators,m_diag,m_spacing,split_ratio, iteration);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        std::cerr << "Guided split in " << elapsed << "[us]" << std::endl;
        return clusters_count;
    }
    Pointset get_point_cloud_clustered()
    {
        Pointset pointset;
        std::vector<Vector> colors;
        for(int k = 0 ; k < m_generators.size();k++)
        {
            double r = (double) rand() / (RAND_MAX);
            double g = (double) rand() / (RAND_MAX);
            double b = (double) rand() / (RAND_MAX);
            colors.push_back(Vector(r,g,b));
        }
        for(int i = 0; i < m_points.size();i++)
        {
            pointset.insert(m_points[i].first,colors[m_vlabels[i]]);
        }
        return pointset;
    }
    const Polyhedron& get_reconstructed()
    {
        return m_triangle_fit.get_mesh();
    }
    void write_csv()
    {
        m_cluster->write_csv();
    }
    // reconstruction 
    void reconstruction(double dist_ratio, double fitting, double coverage, double complexity)
    {
        create_adjacent_edges();
        create_candidate_facets();
        mlp_reconstruction(dist_ratio, fitting, coverage, complexity);
    }
    void create_adjacent_edges()
    {
        if(m_generators.size() == 0)
        {
            if(m_verbose_level > 1)
                std::cout << "No available pole!" << std::endl;
            return;
        }

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        std::vector<Point> dual_points;
        for(int i = 0; i < m_generators.size(); i++)
        {
            Point center;
            // todoquestion
            //if(m_generators_count.size() == m_generators.size())
            //center = pointset_.point(m_generators[i]);
            /*else*/
            center = compute_optimal_point(m_generators_qem[i], pointset_.point(m_generators[i]));

            dual_points.push_back(center);
        }

        IntPairSet adjacent_pairs;
        for(int i = 0; i < pointset_.size(); i++)
        {
            if(m_vlabels.find(i) == m_vlabels.end())
                continue;

            int center_label = m_vlabels[i];
            K_neighbor_search search(m_tree, pointset_.point(i), m_num_knn);

            for(typename K_neighbor_search::iterator it = search.begin(); it != search.end(); it++)
            {
                int nb_index = (it->first).second;

                if(m_vlabels.find(nb_index) != m_vlabels.end())
                {
                    int nb_label = m_vlabels[nb_index];
                    if(center_label != nb_label) 
                        adjacent_pairs.insert(std::make_pair(center_label, nb_label));
                }              
            }
        }
        m_triangle_fit.initialize_adjacent_graph(dual_points, adjacent_pairs,m_bbox,m_diag);

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        if(m_verbose_level > 1)
            std::cerr << "Candidate edge in " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[us]" << std::endl;
    }
    // helper ?
        Point compute_optimal_point(QEM_metric& cluster_qem, Point& cluster_pole)
    {
        // solve Qx = b
        Eigen::MatrixXd qem_mat = cluster_qem.get_4x4_svd_matrix();
        Eigen::VectorXd qem_vec = qem_mat.row(3); // 0., 0., 0., 1.
        Eigen::VectorXd optim(4);

        // check rank
        Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(qem_mat);
        lu_decomp.setThreshold(1e-5);

        // full rank -> direct inverse
        if(lu_decomp.isInvertible())
        {
            optim = lu_decomp.inverse() * qem_vec;
        }
        else
        {   // low rank -> svd pseudo-inverse
            Eigen::JacobiSVD<Eigen::MatrixXd> svd_decomp(qem_mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
            svd_decomp.setThreshold(1e-5);

            optim(0) = cluster_pole.x();
            optim(1) = cluster_pole.y();
            optim(2) = cluster_pole.z();
            optim(3) = 1.;

            optim = optim + svd_decomp.solve(qem_vec - qem_mat * optim);
        }

        Point optim_point(optim(0), optim(1), optim(2));

        return optim_point;
    }
    void update_adjacent_edges(std::vector<float>& adjacent_edges)
    {
        m_triangle_fit.update_adjacent_edges(adjacent_edges);
    }
    void create_candidate_facets()
    {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        m_triangle_fit.create_candidate_facets(); 
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        if(m_verbose_level > 1)
            std::cerr << "Candidate facet in " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[us]" << std::endl;
    }
    void update_candidate_facets(std::vector<float>& candidate_facets, std::vector<float>& candidate_normals)
    {
        m_triangle_fit.update_candidate_facets(candidate_facets, candidate_normals); 
    }

    void mlp_reconstruction(double dist_ratio, double fitting, double coverage, double complexity)
    {
        std::vector<Point> input_point_set;

        std::transform( m_points.begin(), 
                        m_points.end(), 
                        std::back_inserter(input_point_set), 
                        [](const Point_with_index& p) { return p.first; }); 

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        m_triangle_fit.reconstruct(input_point_set, m_spacing, dist_ratio, fitting, coverage, complexity); 
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cerr << "MIP solver in " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[us]" << std::endl;
        auto valid = m_triangle_fit.get_mesh().is_valid();
        if(!valid)
        {
            std::cout<<"Manifold Reconstruction failed, trying with Nonmanifold Reconstruction\n";
            m_triangle_fit.nonmanifold_reconstruct(input_point_set, m_spacing, dist_ratio, fitting, coverage, complexity); 
        }
    }
    void update_fit_surface(std::vector<float>& fit_facets, std::vector<float>& fit_normals)
    {
        m_triangle_fit.update_fit_surface(fit_facets, fit_normals);
    }

    void update_fit_soup(std::vector<float>& fit_soup_facets, std::vector<float>& fit_soup_normals)
    {
        m_triangle_fit.update_fit_soup(fit_soup_facets, fit_soup_normals);
    }
};
}