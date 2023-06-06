#include "types.h"

#include "candidate.h"
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

typedef CGAL::Bbox_3   Bbox;
namespace CGAL {
    class Variational_shape_reconstruction
    {
        private:
        std::vector<QEM_metric>         m_pqems;
        std::vector<QEM_metric>         m_vqems;
        std::vector<std::pair<Point, int>> m_points;
        std::vector<int> m_poles;
        std::vector<Vector> m_normals;
        KNNTree         m_tree;
        int m_num_knn = 7;
        size_t generator_count;

        //RG
            //std::map<int, std::vector<int>>   m_vlabels;
            std::map<int, int>   m_vlabels;
            std::vector<QEM_metric>         m_poles_qem;
            std::vector<int>         m_poles_count;
        //init
        Bbox            m_bbox;
        double          m_diag;
        double          m_spacing;
        public:
        void initialize(Pointset& pointset)
        {		 
            // initialize_point_qem_map
            load_points(pointset);
            initialize_vertex_qem();
            compute_bounding_box();   

            // init kdtree
            m_tree.clear();
            m_tree.insert(m_points.begin(), m_points.end());                                                                       
            // compute average spacing
            //m_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(m_points, 6, CGAL::parameters::point_map(Point_map_pwi()));
            m_spacing=0.01;
        }
        void load_points(Pointset& pointset)
        {
            int num_nb = std::max(6, m_num_knn);
            // note: no const for distance
            for( Pointset::iterator pointset_it = pointset.begin(); pointset_it != pointset.end(); ++ pointset_it )
            {
                auto point = pointset.point(*pointset_it);
                auto normal = pointset.normal(*pointset_it);

                m_points.push_back(std::make_pair(point, std::distance(pointset.begin(),pointset_it)));
                
                K_neighbor_search search(m_tree, point, num_nb);
                KNNDistance tr_dist;

                double avg_dist = 0.;
                for(typename K_neighbor_search::iterator it = search.begin(); it != search.end(); it++)
                    avg_dist += tr_dist.inverse_of_transformed_distance(it->second);

                avg_dist = avg_dist / (double)num_nb;
                
                QEM_metric pqem = compute_qem_for_point(point, normal, avg_dist * avg_dist);
                m_pqems.push_back(pqem);
        
            }
            std::cout << "Number of points: " << m_points.size() << std::endl;
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
            boost::function<Point(std::pair<Point, int>&)> pwi_it_to_point_it = boost::bind(&std::pair<Point, int>::first, _1);
            m_bbox = CGAL::bbox_3(  boost::make_transform_iterator(m_points.begin(), pwi_it_to_point_it), 
                                    boost::make_transform_iterator(m_points.end(), pwi_it_to_point_it));
            m_diag = std::sqrt(CGAL::squared_distance(Point(m_bbox.min(0), m_bbox.min(1), m_bbox.min(2)),
                                                    Point(m_bbox.max(0), m_bbox.max(1), m_bbox.max(2))));
            std::cout << "Diagonal of bounding box: " << m_diag << std::endl;
        }


        void init_random_poles(int num_poles)
        {
            generator_count = num_poles;
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

            std::vector<int> num_range(m_points.size());
            std::iota(num_range.begin(), num_range.end(), 0);

            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(num_range.begin(), num_range.end(), g);

            std::set<int> selected_indices;
            for(int i = 0; i < num_poles; i++)
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

            std::cout << "m_vlabels size: " << m_vlabels.size() << std::endl;
            std::cout << "m_points size: " << m_points.size() << std::endl;
            std::cout<< "iterationrs "<<k<<"\n";
            std::cout<< "iterationrs2 "<<p<<"\n";
            assert(m_vlabels.size() == m_points.size());
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cerr << "\nRegion growing in " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000 << "[ms]" << std::endl;
            savePs(m_points,m_vlabels, generator_count,"output.txt");
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
                center = m_points[m_poles[i]].first;
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
    };
} // namespace CGAL