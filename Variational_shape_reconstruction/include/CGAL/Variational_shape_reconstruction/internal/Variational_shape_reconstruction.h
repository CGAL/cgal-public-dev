#include "types.h"

#include "candidate.h"
#include "trianglefit.h"
#include "pqueue.h"
#include "io.h"
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Splitters.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include "metrics.h"
#include <CGAL/bounding_box.h>
#include <CGAL/compute_average_spacing.h>

#include <Eigen/Dense>
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/iterator/transform_iterator.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <random>



typedef std::pair<Point, std::size_t>                                               Point_with_index;
typedef std::vector<Point_with_index>                                               PwiList;
typedef CGAL::First_of_pair_property_map<Point_with_index>                          Point_map_pwi;
// knntree
typedef CGAL::Search_traits_3<Kernel>                                               Traits_base;
typedef CGAL::Search_traits_adapter<Point_with_index,Point_map_pwi, Traits_base>    Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>                                  K_neighbor_search;
typedef typename K_neighbor_search::Tree                                            KNNTree;
typedef typename K_neighbor_search::Distance                                        KNNDistance;

//Pqueue
typedef qem::Candidate<int>                                          CCandidate;
typedef qem::Candidate_more<CCandidate>                              More;
typedef qem::Custom_priority_queue<CCandidate, More>                 PQueue;
namespace qem {
typedef typename std::unordered_set<std::pair<int, int>, HashPairIndex, EqualPairIndex>     IntPairSet; 

typedef CGAL::Bbox_3   Bbox;

    class Variational_shape_reconstruction
    {
        private:
            //qem
            std::vector<QEM_metric> m_pqems;
            std::vector<QEM_metric> m_vqems;

            //generators
            size_t m_generator_count;
            std::vector<int> m_poles;

            // geometry
            std::vector<std::pair<Point, size_t>> m_points;
            std::vector<Vector> m_normals;

            // knntree
            KNNTree m_tree;
            int m_num_knn = 7;
            

            //Region growing
            std::map<int, int> m_vlabels;
            std::vector<QEM_metric> m_poles_qem;
            std::vector<int> m_poles_count;

            //init
            Bbox            m_bbox;
            double          m_diag;
            double          m_spacing;

            TriangleFit              m_triangle_fit;


        public:
        Variational_shape_reconstruction(int generator_count) : m_generator_count(generator_count) {}
        void initialize(Pointset& pointset)
        {		 
            load_points(pointset);
            compute_bounding_box();   

            // init kdtree
            m_tree.clear();
            m_tree.insert(m_points.begin(), m_points.end());     
                                                              
            // compute average spacing
            m_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(m_points, 6, CGAL::parameters::point_map(Point_map_pwi()));            
            initialize_qem_map();
            initialize_vertex_qem();    
            init_random_poles();
        }
        void load_points(Pointset& pointset)
        {
            // note: no const for distance
            for( Pointset::iterator pointset_it = pointset.begin(); pointset_it != pointset.end(); ++ pointset_it )
            {
                const auto point = pointset.point(*pointset_it);
                const auto normal = pointset.normal(*pointset_it);

                m_points.push_back(std::make_pair(point, std::distance(pointset.begin(),pointset_it)));
                m_normals.push_back(normal);
        
            }
            std::cout << "Number of points: " << m_points.size() << std::endl;
        }
        void initialize_qem_map()
        {
             int num_nb = std::max(6, m_num_knn);
            for(int i = 0; i < m_points.size(); i++)
            {
                K_neighbor_search search(m_tree, m_points[i].first, num_nb);
                KNNDistance tr_dist;

                double avg_dist = 0.;
                for(typename K_neighbor_search::iterator it = search.begin(); it != search.end(); it++)
                    avg_dist += tr_dist.inverse_of_transformed_distance(it->second);

                avg_dist = avg_dist / (double)num_nb;
                
                QEM_metric pqem = compute_qem_for_point(m_points[i].first, m_normals[i], avg_dist * avg_dist);
                m_pqems.push_back(pqem);
            }
        }
        void initialize_vertex_qem()
        {
            m_vqems.clear();
            for(int i = 0; i < m_points.size(); i++)
            {
                K_neighbor_search search(m_tree, m_points[i].first, m_num_knn);
                QEM_metric vqem;

                for(typename K_neighbor_search::iterator it = search.begin(); it != search.end(); it++)
                    vqem = vqem + m_pqems[(it->first).second];

                m_vqems.push_back(vqem);
            }
        }
        void compute_bounding_box()
        {
            // find bounding box
            boost::function<Point(std::pair<Point, size_t>&)> pwi_it_to_point_it = boost::bind(&std::pair<Point, size_t>::first, _1);
            m_bbox = CGAL::bbox_3(  boost::make_transform_iterator(m_points.begin(), pwi_it_to_point_it), 
                                    boost::make_transform_iterator(m_points.end(), pwi_it_to_point_it));
            m_diag = std::sqrt(CGAL::squared_distance(Point(m_bbox.min(0), m_bbox.min(1), m_bbox.min(2)),
                                                    Point(m_bbox.max(0), m_bbox.max(1), m_bbox.max(2))));
            std::cout << "Diagonal of bounding box: " << m_diag << std::endl;
        }


        void init_random_poles()
        {
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

            std::vector<int> num_range(m_points.size());
            std::iota(num_range.begin(), num_range.end(), 0);

            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(num_range.begin(), num_range.end(), g);

            std::set<int> selected_indices;
            for(int i = 0; i < m_generator_count; i++)
                selected_indices.insert(num_range[i]);
                
            for(auto &elem: selected_indices) {
                m_poles.push_back(elem);
            }   

            std::cout << "Number of random poles: " << m_poles.size() << std::endl;

            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cerr << "Random poles in " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[us]" << std::endl;
            
        }

        void region_growing(bool flag_dist)
        {
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            m_vlabels.clear();
            m_poles_qem.clear();

            PQueue growing_queue;
            
            // init seed triangles
            for(int label = 0; label < m_poles.size(); label++)
            {
                int index = m_poles[label];
                m_vlabels[index] =  label;
                m_poles_qem.push_back(m_vqems[index]);
                add_candidates(growing_queue,index,label, flag_dist);
            }
            int k =0;
            int p =0;
            while(!growing_queue.empty())
            {
                
                const CCandidate candidate = growing_queue.top();
                growing_queue.pop();
                const int index = candidate.handle();
                const int label = candidate.index();
                if(m_vlabels.find(index) != m_vlabels.end())
                {
                    p++;
                    continue;
                }
                k++;
                m_vlabels[index] = label;
                m_poles_qem[label] = m_poles_qem[label] + m_vqems[index]; 
                add_candidates(growing_queue,index,label, flag_dist);
            }

            //std::cout << "m_vlabels size: " << m_vlabels.size() << std::endl;
            //std::cout << "m_points size: " << m_points.size() << std::endl;
            //std::cout<< "iterations "<<k<<"\n";
            //std::cout<< "iterations2 "<<p<<"\n";
            assert(m_vlabels.size() == m_points.size());
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cerr << "\nRegion growing in " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000 << "[ms]" << std::endl;
            savePs(m_points,m_vlabels,m_poles.size(),"output.txt");
        }
        void add_candidates(PQueue &growing_queue,const int index,const int label,const bool flag_dist)
        {
            
            K_neighbor_search search(m_tree, m_points[index].first, m_num_knn);
            for(typename K_neighbor_search::iterator it = search.begin(); it != search.end(); it++)
            {
                const int nb_index = (it->first).second;

                if(m_vlabels.find(nb_index) == m_vlabels.end() )
                {
                    const double loss = compute_collapse_loss(nb_index, label, flag_dist);
                    growing_queue.push(CCandidate(nb_index, label, loss));
                }
            }
        }

        bool update_poles()
        {
            std::cout << "Updating poles..." << std::endl;

            std::vector<Point>    optimal_points;
            std::vector<double> dists;
            std::vector<int> old_poles;

            for(int i = 0; i < m_poles.size(); i++)
            {
                Point center;
                // is it for multicover ?
                /*if(m_poles_count.size() == m_poles.size() && m_poles_count[i] == 1)
                center = m_points[m_poles[i]].first;
                else*/
                center = compute_optimal_point(m_poles_qem[i], m_points[m_poles[i]].first);
                optimal_points.push_back(center);
                dists.push_back(1e20);
                old_poles.push_back(m_poles[i]);
            }
            
            for(int i = 0; i < m_points.size(); i++) 
            {
                if(m_vlabels.find(i) == m_vlabels.end())
                    continue;

                int label = m_vlabels[i];
                double dist = CGAL::squared_distance(optimal_points[label], m_points[i].first);

                if(dist < dists[label])
                {
                    m_poles[label] = i;
                    dists[label] = dist;
                }
            }

            // chech change
            for(int i = 0; i < m_poles.size(); i++)
            {
                if(m_poles[i] != old_poles[i])
                    return true;
            }

            std::cout << "Region growing converges!" << std::endl;
            return false;
        }
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
    size_t guided_split_clusters(double split_ratio, size_t iteration) // batch splitting
    {
        std::cout << "Begin guided split..." << std::endl;

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        double split_thresh = m_diag * split_ratio;
        split_thresh = split_thresh * split_thresh; // square distance
        split_thresh = split_thresh * m_num_knn * m_spacing * m_spacing;
        std::cout << "  split threshold: " << split_thresh << std::endl;

        // compute error for each center
        std::vector<double> qem_errors(m_poles.size(), 0.);
        std::vector<int> vert_indices(m_poles.size(), -1);

        for(int i = 0; i < m_points.size(); i++) 
        {
            if(m_vlabels.find(i) == m_vlabels.end())
                continue;

            int center_ind = m_vlabels[i];
            double error = compute_minimum_qem_error(m_points[m_poles[center_ind]].first, m_vqems[i]);

            if(error > qem_errors[center_ind])
            {
                qem_errors[center_ind] = error;
                vert_indices[center_ind] = i;
            }
        }

        // split centers exceeding max error
        std::vector<int> new_poles;

        for(int i = 0; i < m_poles.size(); i++)
        {

            if(qem_errors[i] > split_thresh && vert_indices[i] != m_poles[i])
                new_poles.push_back(vert_indices[i]);
        }

    std::cout << "Found " << new_poles.size() << " new poles!" << std::endl;

        // merge close poles
        std::set<int, std::greater<int> > duplicate_poles;
        double dist_thresh = m_spacing * 5.;
        dist_thresh = dist_thresh * dist_thresh;

        for(int i = 0; i < new_poles.size(); i++)
        {
            for(int j = i + 1; j < new_poles.size(); j++)
            {
                if(duplicate_poles.find(j) == duplicate_poles.end() && 
                   CGAL::squared_distance(m_points[new_poles[i]].first, m_points[new_poles[j]].first) < dist_thresh)
                {
                    duplicate_poles.insert(j);
                }
            }
        }

        for(auto& elem: duplicate_poles)
        {
            new_poles.erase(new_poles.begin() + elem);
        }

        std::cout << "Remove " << duplicate_poles.size() << " duplicated poles!" << std::endl;

        // insert new poles
        for(int i = 0; i < new_poles.size(); i++)
        {
            m_poles.push_back(new_poles[i]);

        }

        std::cout << "Finally found " << new_poles.size() << " new poles!" << std::endl;

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	    std::cerr << "Guided split in " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[us]" << std::endl;
        savePs(m_points,m_vlabels,m_poles.size(),"output_guided_"+std::to_string(iteration)+".txt");
        return new_poles.size();
    }
        double compute_collapse_loss(const int index,const int label,const bool flag) //eq 6
        {
            const double m_qem_weight=1.;
            const double m_dist_weight=1.;
            
            const double qem_cost = compute_minimum_qem_error(m_points[m_poles[label]].first, m_vqems[index]);

            const double cost = m_qem_weight * qem_cost;
            // if flag
            const double dist_cost = m_num_knn * CGAL::squared_distance(m_points[m_poles[label]].first, m_points[index].first);       

            return cost + m_dist_weight * dist_cost;
        }
        double kernel_function(double dist, double h)
        {
            double value = 0.;
            double ratio = dist / h;

            if(ratio >= 1.)
                return value;

            value = std::pow(1. - std::pow(ratio, 2), 4);

            return value;
        } 
        double compute_minimum_qem_error(Point center_point, QEM_metric& query_qem)
        {
            Eigen::VectorXd center_vec(4);
            center_vec << center_point.x(), center_point.y(), center_point.z(), 1.;

            double cost = center_vec.transpose() * query_qem.get_4x4_matrix() * center_vec;

            return std::abs(cost);
        }
        QEM_metric compute_qem_for_point(const Point& query,const Vector& normal,const double &area)
        {
            QEM_metric qem;
            qem.init_qem_metrics_face(area, query, normal);
            return qem;
        }
            // adjacente
       void create_adjacent_edges()
        {
            if(m_poles.size() == 0)
            {
                std::cout << "No available pole!" << std::endl;
                return;
            }

            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

            std::vector<Point> dual_points;
            for(int i = 0; i < m_poles.size(); i++)
            {
                Point center;

                if(m_poles_count.size() == m_poles.size() && m_poles_count[i] == 1)
                    center = m_points[m_poles[i]].first;
                else
                    center = compute_optimal_point(m_poles_qem[i], m_points[m_poles[i]].first);

                dual_points.push_back(center);
            }

            IntPairSet adjacent_pairs;
            for(int i = 0; i < m_points.size(); i++)
            {
                if(m_vlabels.find(i) == m_vlabels.end())
                    continue;

                int center_label = m_vlabels[i];
                K_neighbor_search search(m_tree, m_points[i].first, m_num_knn);

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

            m_triangle_fit.initialize_adjacent_graph(dual_points, m_poles_qem, adjacent_pairs);

            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cerr << "Candidate edge in " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[us]" << std::endl;
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
        }
        void update_fit_surface(std::vector<float>& fit_facets, std::vector<float>& fit_normals)
        {
            m_triangle_fit.update_fit_surface(fit_facets, fit_normals);
        }

        void update_fit_soup(std::vector<float>& fit_soup_facets, std::vector<float>& fit_soup_normals)
        {
            m_triangle_fit.update_fit_soup(fit_soup_facets, fit_soup_normals);
            m_triangle_fit.save_trianglefit_mesh("filename.off");
        }
    };
}