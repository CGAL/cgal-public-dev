#include "types.h"
//qem
#include "metrics.h"
#include "candidate.h"
#include "pqueue.h"

#include <CGAL/bounding_box.h>
#include <CGAL/compute_average_spacing.h>
// knn
#include <CGAL/Orthogonal_k_neighbor_search.h>

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/iterator/transform_iterator.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <random>

#include "helper_metrics.h"
#include "io.h"

typedef std::pair<Point, std::size_t>                                               Point_with_index;
typedef std::vector<Point_with_index>                                               PwiList;
typedef CGAL::First_of_pair_property_map<Point_with_index>                          Point_map_pwi;

// knntree
typedef CGAL::Search_traits_3<Kernel>                                               Traits_base;
typedef CGAL::Search_traits_adapter<Point_with_index, Point_map_pwi, Traits_base>   Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>                                  K_neighbor_search;
typedef typename K_neighbor_search::Tree                                            KNNTree;
typedef typename K_neighbor_search::Distance                                        KNNDistance;

//Pqueue
typedef qem::Candidate<int>                                          CCandidate;
typedef qem::Candidate_more<CCandidate>                              More;
typedef qem::Custom_priority_queue<CCandidate, More>                 PQueue;

namespace qem
{   
    class Clustering
    {   
        public:
        Clustering()
        {
        }

        Clustering(const Pointset& pointset,
        size_t num_knn,
        double euclidean_distance_weight,
        VERBOSE_LEVEL verbose_level)
        {
            m_point_set = pointset;
            m_num_knn = num_knn;
            m_dist_weight = euclidean_distance_weight;
            m_verbose_level = verbose_level;
            csv_writer =std::make_shared<DataWriter>(pointset.size());
        }

        /// @brief Compute the qem for each point based on the k nearest neighbor neighbors
        // fixme: rather based on the normal and average distance to neighbors!
        void initialize_qem_per_point(const KNNTree& m_tree)
        {
            // init vector of qems
            m_pqems.clear();

            int num_nb = std::max(6, m_num_knn); // fixme

            for(Pointset::const_iterator it = m_point_set.begin(); it != m_point_set.end(); it++)
            {
                auto point = m_point_set.point(*it);
                K_neighbor_search search(m_tree, point, num_nb);
                KNNDistance tr_dist;

                double avg_dist = 0.;
                for(typename K_neighbor_search::iterator it = search.begin(); it != search.end(); it++)
                    avg_dist += tr_dist.inverse_of_transformed_distance(it->second);

                avg_dist = avg_dist / (double)num_nb;
                
                QEM_metric pqem = compute_qem_for_point(point, m_point_set.normal(*it), avg_dist * avg_dist);
                m_pqems.push_back(pqem);
            }
        }

        /// @brief Compute the connected components
        void compute_connected()
        {
            size_t point_count = m_graph.size(); 
            m_visited.resize(point_count);
            m_component.resize(point_count);
            std::fill(m_visited.begin(), m_visited.end(), false);
            std::fill(m_component.begin(), m_component.end(), 0);
            
            for(int i = 0 ; i < point_count;i++)
            {
                if(!m_visited[i])
                {
                    visit(i, component_count);
                    component_count++;
                }
            }
            
            if(m_verbose_level == VERBOSE_LEVEL::HIGH)
            {
                std::ofstream clustering_connected;
                clustering_connected.open("clustering_connected.ply");

                clustering_connected << "ply\n"
                            << "format ascii 1.0\n"
                            << "element vertex " << m_point_set.size() << "\n"
                            << "property float x\n"
                            << "property float y\n"
                            << "property float z\n"
                            << "property uchar red\n"
                            << "property uchar green\n"
                            << "property uchar blue\n"
                            << "end_header\n";
                std::vector<Vector> colors;
                for(int i = 0 ; i < component_count; i++)
                {
                    double r = (double) rand() / (RAND_MAX);
                    double g = (double) rand() / (RAND_MAX);
                    double b = (double) rand() / (RAND_MAX);
                    colors.push_back(Vector(r,g,b));
                }
                int point_index =0;
                for(Pointset::const_iterator it = m_point_set.begin(); it != m_point_set.end(); ++ it)
                {

                    auto point = m_point_set.point(*it);
                    clustering_connected << point.x() << " " << point.y() << " " << point.z() << " ";
                    auto normal = colors[m_component[point_index]];
                    clustering_connected << static_cast<int>(255*normal.x()) << " " << static_cast<int>(255*normal.y()) << " " << static_cast<int>(255*normal.z()) << std::endl;

                    point_index++;
                }
                clustering_connected.close();
            }
            if(m_verbose_level != VERBOSE_LEVEL::LOW)
            {
                std::cout<< "m_component.size() " << m_component.size() << std::endl;
                std::cout<< "point_count " << point_count << std::endl;
            }

        }
        /// @brief visit function that explore the graph of neighbors using a 
        /// queue to assign a component index to each point
        /// @param current_idx index of the current point
        /// @param component_idx index of the current connected component
        void visit(const int current_idx, const int component_idx)
        {
            std::vector<int> queue;
            queue.push_back(current_idx);
            
            while(!queue.empty())
            {
                int idx = queue.back();
                queue.pop_back();
                m_component[idx] = component_idx;
                m_visited[idx] = true;
                for(auto neighbors_idx : m_graph[idx])
                {
                    if(!m_visited[neighbors_idx])
                    {
                        queue.push_back(neighbors_idx);
                    }
                }
            }
        }
        /// @brief Compute the qem for a point weighted by the area of the face
        /// @param query the point to compute the qem for
        /// @param normal of the face
        /// @param area of the face
        /// @return the qem computed
        QEM_metric compute_qem_for_point(const Point& query, const Vector& normal, const double &area)
        {
            QEM_metric qem;
            qem.init_qem_metrics_face(area, query, normal);
            return qem;
        }
        /// @brief Compute the sum of the qem neighbor points  for each point in m_vqems
        /// Also build the graph of neighbors 
        /// @param m_tree the knn tree
        /// @param generators the list of indices of generators
        /// // fixme: the connected components and generator addition should be performed in a separate function
        void initialize_qem_per_vertex(const KNNTree& m_tree, std::vector<int>& generators)
        {
            m_vqems.clear();
            int point_index = 0;
            for(Pointset::const_iterator it = m_point_set.begin(); it != m_point_set.end(); it++, point_index++)
            {
                K_neighbor_search search(m_tree, m_point_set.point(*it), m_num_knn);

                std::vector<int> neighbors;

                // init with qem of point itself 
                QEM_metric vqem = m_pqems[point_index];

                for(typename K_neighbor_search::iterator it = search.begin(); it != search.end(); it++)
                {
                    auto neighbor_idx = (it->first).second;
                    vqem = vqem + m_pqems[neighbor_idx];
                    neighbors.push_back(neighbor_idx);
                }
                m_graph.push_back(neighbors);

                m_vqems.push_back(vqem);
            }

            // compute connected components
            compute_connected();

            if(m_verbose_level != VERBOSE_LEVEL::LOW)
            {
                std::cout << "#connected components: " << component_count - 1 << std::endl;
                std::cout << "#generators: " << generators.size() << std::endl;
            }
            for(int i = 0 ; i < component_count;i++)
            {
                bool found = false;
                for(int j = 0 ; j < generators.size();j++)
                {
                    if(i == m_component[generators[j]])
                    {
                        found = true;
                    }
                }
                if(!found)
                {
                    // Add a new generator in the connected component without generator
                    auto it = std::find(m_component.begin(), m_component.end(), i);
                    int index = it - m_component.begin();
                    generators.push_back(index);
                    m_component[index]=i;
                    if(m_verbose_level != VERBOSE_LEVEL::LOW)
                        std::cout<<"added generator > "<<index<<" connected component "<< i<<std::endl;
                    
                }
            }

            if(m_verbose_level == VERBOSE_LEVEL::HIGH)
            {
                // Create a pointcloud of the graph of neighbors so that
                // each point is connected to each of his neighbors
                std::ofstream neighbors_graph;
                neighbors_graph.open("neighbors_graph.ply");

                std::size_t sum = 0;
                for (auto &&i : m_graph) {
                    sum += i.size();
                }

                neighbors_graph << "ply\n"
                        << "format ascii 1.0\n"
                        << "element vertex " << m_point_set.size() << "\n"
                        << "property float x\n"
                        << "property float y\n"
                        << "property float z\n"
                        << "element face " << sum << "\n"
                        << "property list uchar int vertex_indices\n"
                        << "end_header\n";

                for(Pointset::const_iterator it = m_point_set.begin(); it != m_point_set.end(); ++ it)
                {
                    auto point = m_point_set.point(*it);
                    neighbors_graph << point.x() << " " << point.y() << " " << point.z() << std::endl;
                }
                
                for(int i = 0; i < m_graph.size(); i++)
                {
                    for(int j = 0; j < m_graph[i].size(); j++)
                    {
                        neighbors_graph << "2 "<<i<<" "<<m_graph[i][j]<<"\n";
                    }
                }
                neighbors_graph.close();
            }

        }
        int get_component_count()
        {
            return component_count;
        }
        std::vector<int> get_generators_per_component(std::vector<int>& generators)
        {
            std::vector<int> generators_per_component(component_count,0);
            for(int j = 0 ; j < generators.size();j++)
            {
                generators_per_component[m_component[generators[j]]]++;
            }

            return generators_per_component;
        }

        /// @brief Partition: find the best generator for each point, and update cluster QEM
        /// via region growing and the cost function compute_collapse_loss
        /// @param m_vlabels 
        /// @param generators_qem 
        /// @param generators 
        /// @param flag_dist 
        void partition(std::map<int, int>& m_vlabels,
        std::vector<QEM_metric>& generators_qem,
        const std::vector<int>& generators,
        const bool flag_dist)
        {
            PQueue growing_queue;
            
            // init seed triangles
            for(int label = 0; label < generators.size(); label++)
            {
                int index = generators[label];
                m_vlabels[index] =  label;
                generators_qem.push_back(m_vqems[index]);
                add_candidates(growing_queue, index, label, flag_dist, m_vlabels, generators);
            }

            while(!growing_queue.empty())
            {
                const CCandidate candidate = growing_queue.top();
                growing_queue.pop();
                const int index = candidate.handle();
                const int label = candidate.index();

                if(m_vlabels.find(index) != m_vlabels.end())
                    continue;

                m_vlabels[index] = label;
                generators_qem[label] = generators_qem[label] + m_vqems[index]; 
                add_candidates(growing_queue, index, label, flag_dist, m_vlabels, generators);
            }
        }
        /// @brief Add the generators candidates to the priority queue
        /// @param growing_queue 
        /// @param index index of the current point
        /// @param label index of the generator associated to the point
        /// @param flag_dist 
        /// @param m_vlabels 
        /// @param generators 
        void add_candidates(PQueue &growing_queue,
        const int index,
        const int label,
        const bool flag_dist,
        const std::map<int, int>& m_vlabels,
        const std::vector<int>& generators)
        {
            for(const auto nb_index : m_graph[index])
            {
                if(m_vlabels.find(nb_index) == m_vlabels.end() )
                {
                    const double loss = compute_collapse_loss(nb_index, label, flag_dist, generators);
                    growing_queue.push(CCandidate(nb_index, label, loss));
                }
            }
        }

        /// @brief Compute the collapse loss using the qem cost and weighted by the euclidean distance
        /// @param index index of the current point
        /// @param label index of the generator associated to the point
        /// @param flag flag to use the euclidean distance
        /// @param generators 
        /// @return the cost
        double compute_collapse_loss(const int index, const int label, const bool flag, const std::vector<int>& generators) //eq 6
        {
            const double qem_cost = compute_minimum_qem_error(m_point_set.point(generators[label]), m_vqems[index]);
            double cost = qem_cost; 

            if(flag)
            {
                double dist_cost = m_num_knn * CGAL::squared_distance(m_point_set.point(generators[label]), m_point_set.point(index));
                cost = cost + m_dist_weight * dist_cost;
            }
            return cost;
        }

        /// @brief Compute the minimun qem error from a center point with a query_qem
        /// @param center_point 
        /// @param query_qem 
        /// @return the minimum qem error
        /// FIXME: rename - this is just the QEM error for a give point
        double compute_minimum_qem_error(const Point& center_point, QEM_metric& query_qem)
        {
            Eigen::VectorXd center_vec(4);
            center_vec << center_point.x(), center_point.y(), center_point.z(), 1.;

            double error = center_vec.transpose() * query_qem.get_4x4_matrix() * center_vec;

            return std::abs(error); // fixme: check that it is not negative instead
        }

        /// @brief Update the generators 
        /// @param m_vlabels 
        /// @param generators_qem 
        /// @param generators 
        /// @return a Boolean : true if generators changed and false otherwise 
        bool update_generators(std::map<int, int>& m_vlabels,
        std::vector<QEM_metric>& generators_qem,
        std::vector<int>& generators)
        {
            if(m_verbose_level != VERBOSE_LEVEL::LOW)
                std::cout << "Update generators" << std::endl;

            std::vector<Point> optimal_points;
            std::vector<double> dists;
            std::vector<int> old_generators;

            for(int i = 0; i < generators.size(); i++)
            {
                Point center = compute_optimal_point(generators_qem[i], m_point_set.point(generators[i]));
                optimal_points.push_back(center);
                dists.push_back(1e20); // fixme
                old_generators.push_back(generators[i]);
            }
            
            for(int i = 0; i < m_point_set.size(); i++) 
            {
                if(m_vlabels.find(i) == m_vlabels.end())
                    continue;

                int label = m_vlabels[i];
                double dist = CGAL::squared_distance(optimal_points[label], m_point_set.point(i));

                if(dist < dists[label])
                {
                    generators[label] = i;
                    dists[label] = dist;
                }
            }

            // compte errors
            std::vector<double> qem_errors(generators.size(), 0.);
            for(int i = 0; i < m_point_set.size(); i++) 
            {
                if(m_vlabels.find(i) == m_vlabels.end())
                    continue;

                int center_ind = m_vlabels[i];

                // fixme: it should be the QEM error of the optimal generator point! not of the input point coinciding with the generator
                double error = compute_minimum_qem_error(m_point_set.point(generators[center_ind]), m_vqems[i]); 
                
                
                csv_writer->addErrorPoints(i,error);

                if(error > qem_errors[center_ind])
                {
                    qem_errors[center_ind] = error;
                }
            }

            double error = std::numeric_limits<double>::min();
            for(int i = 0 ; i < qem_errors.size();i++)
            {
                csv_writer->addWorstErrorGenerator(i,qem_errors[i]);
                error = std::max(error, qem_errors[i]);
            }

            auto mean = std::accumulate(qem_errors.begin(), qem_errors.end(), 0.);
            mean /= qem_errors.size();
            csv_writer->addWorstErrorGenerator(error);
            csv_writer->addMeanErrorGenerator(mean);
            csv_writer->setGenerator(generators.size());

            if(m_verbose_level == VERBOSE_LEVEL::HIGH)
            {
                std::ofstream clustering_by_iteration;
                clustering_by_iteration.open("output/clustering_"+std::to_string(m_id++)+".ply");

                clustering_by_iteration << "ply\n"
                            << "format ascii 1.0\n"
                            << "element vertex " << m_point_set.size() << "\n"
                            << "property float x\n"
                            << "property float y\n"
                            << "property float z\n"
                            << "property uchar red\n"
                            << "property uchar green\n"
                            << "property uchar blue\n"
                            << "end_header\n";
                std::vector<Point> points;
                for(Pointset::const_iterator it = m_point_set.begin(); it != m_point_set.end(); ++ it)
                {
                    auto point = m_point_set.point(*it);
                    points.push_back(point);
                    clustering_by_iteration << point.x() << " " << point.y() << " " << point.z() << " ";
                    auto normal = m_point_set.normal(*it);
                    clustering_by_iteration << static_cast<int>(255*normal.x()) << " " << static_cast<int>(255*normal.y()) << " " << static_cast<int>(255*normal.z()) << std::endl;
                }

                clustering_by_iteration.close();
            }

            // check change
            for(int i = 0; i < generators.size(); i++)
            {
                if(generators[i] != old_generators[i])
                    return true;
            }

            if(m_verbose_level != VERBOSE_LEVEL::LOW)
                std::cout << "Region growing converged" << std::endl;
            
            return false;
        }
        /// @brief Compute optimal point using either SVD or the direct inverse
        /// @param cluster_qem 
        /// @param cluster_pole 
        /// @return the optimal point
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
        /// @brief Splits the cluster if the qem error is more than a threshold split_thresh
        /// @param m_vlabels 
        /// @param generators 
        /// @param m_diag diagonal of the aabb box
        /// @param m_spacing average spacing between the points 
        /// @param split_ratio user defined parameter for the split
        /// @param iteration 
        /// @return total of generators added
        size_t guided_split_clusters(
        std::map<int, int>& m_vlabels,
        std::vector<int>& generators,
        const double m_diag,
        const double m_spacing,
        const double split_ratio,
        size_t iteration) // batch splitting
        {

            double split_thresh = m_diag * split_ratio;
            split_thresh = split_thresh * split_thresh; // square distance
            split_thresh = split_thresh * m_num_knn * m_spacing * m_spacing;

            // compute error for each center
            std::vector<double> qem_errors(generators.size(), 0.);
            std::vector<int> vert_indices(generators.size(), -1);

            std::map<int,double> generator_worst_error;
            for(int i = 0; i < m_point_set.size(); i++) 
            {
                if(m_vlabels.find(i) == m_vlabels.end())
                    continue;

                int center_ind = m_vlabels[i];
                double error = compute_minimum_qem_error(m_point_set.point(generators[center_ind]), m_vqems[i]); 

                if(error > qem_errors[center_ind])
                {
                    qem_errors[center_ind] = error;
                    vert_indices[center_ind] = i;
                }
            }

            // split centers exceeding max error
            std::vector<int> new_generators;

            for(int i = 0; i < generators.size(); i++)
            {

                if(qem_errors[i] > split_thresh && vert_indices[i] != generators[i])
                    new_generators.push_back(vert_indices[i]);
            }

            std::cout << "Found " << new_generators.size() << " new generators!" << std::endl;

            // merge close generators
            std::set<int, std::greater<int> > duplicate_generators;
            double dist_thresh = m_spacing * 5.;
            dist_thresh = dist_thresh * dist_thresh;

            for(int i = 0; i < new_generators.size(); i++)
            {
                for(int j = i + 1; j < new_generators.size(); j++)
                {
                    if(duplicate_generators.find(j) == duplicate_generators.end() && 
                        CGAL::squared_distance(m_point_set.point(new_generators[i]), m_point_set.point(new_generators[j])) < dist_thresh)
                    {
                        duplicate_generators.insert(j);
                    }
                }
            }

            for(auto& elem: duplicate_generators)
            {
                new_generators.erase(new_generators.begin() + elem);
            }

            std::cout << "Remove " << duplicate_generators.size() << " duplicated generators!" << std::endl;

            // insert new generators
            for(int i = 0; i < new_generators.size(); i++)
            {
                generators.push_back(new_generators[i]);

            }

            std::cout << "Finally found " << new_generators.size() << " new generators!" << std::endl;


            return new_generators.size();
        }
        /// @brief write some data to csv files and print worst error at each iteration 
        void write_csv()
        {
            csv_writer->writeDataErrorGeneratorsToCSV("error_generators.csv");
            csv_writer->writeDataErrorPointsToCSV("error_points.csv");
            csv_writer->printWorst();
        }
        void set_distance_weight(double dist_weight)
        {
            m_dist_weight = dist_weight;
        }

            private:
                Pointset m_point_set;
                int m_num_knn = 12;
                double m_dist_weight = 0.1;

                // qem
                std::vector<QEM_metric> m_pqems; // qem per point
                std::vector<QEM_metric> m_vqems; // diffused qem per point
                std::vector<std::vector<int> > m_graph; // neighborhood graph (indices)
                std::vector<bool> m_visited;
                std::vector<int> m_component;
                int component_count = 0; // fixme
                VERBOSE_LEVEL m_verbose_level;

                // csv
                int m_id = 0;
                std::shared_ptr<DataWriter> csv_writer;
            
    };
}