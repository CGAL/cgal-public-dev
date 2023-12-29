#include "types.h"

// qem
#include "qem.h"
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


// knntree
typedef std::pair<Point, std::size_t> Point_with_index;
typedef std::vector<Point_with_index> PwiList;
typedef CGAL::First_of_pair_property_map<Point_with_index> Point_map_pwi;
typedef CGAL::Search_traits_3<Kernel>                                               Traits_base;
typedef CGAL::Search_traits_adapter<Point_with_index, Point_map_pwi, Traits_base>   Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>                                  K_neighbor_search;
typedef typename K_neighbor_search::Tree                                            KNNTree;
typedef typename K_neighbor_search::Distance                                        KNNDistance;

//Pqueue
typedef qem::Candidate<int> CCandidate;
typedef qem::Candidate_more<CCandidate> More;
typedef qem::Custom_priority_queue<CCandidate, More> PQueue;

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
            csv_writer = std::make_shared<DataWriter>(pointset.size());
        }

        /// @brief Compute the qem for each point based on the k nearest neighbor neighbors
        // fixme: rather based on the normal and average distance to neighbors!
        // TODO: add function to estimate normals
        void initialize_qem_per_point(const KNNTree& tree)
        {
            // init vector of qems
            m_pqems.clear();

            int num_nb = std::max(6, m_num_knn); // fixme

            for(Pointset::const_iterator it = m_point_set.begin(); it != m_point_set.end(); it++)
            {
                auto point = m_point_set.point(*it);
                K_neighbor_search search(tree, point, num_nb);
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
        void initialize_qem_per_vertex(const KNNTree& m_tree) // std::vector<int>& generators
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

            /*

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
            */

        }
        /*
        int get_component_count()
        {
            return component_count;
        }
        */

        /*
        std::vector<int> get_generators_per_component(std::vector<int>& generators)
        {
            std::vector<int> generators_per_component(component_count,0);
            for(int j = 0 ; j < generators.size();j++)
            {
                generators_per_component[m_component[generators[j]]]++;
            }

            return generators_per_component;
        }
        */

        /// @brief Partition: find the best generator for each point, and update cluster QEM
        /// via region growing and the cost function compute_collapse_loss
        /// @param m_vlabels 
        /// @param generators_qem 
        /// @param generators 
        /// @param flag_dist 
        void partition(std::map<int, int>& m_vlabels,
            std::vector<Generator>& generators,
            const bool flag_dist)
        {
            PQueue pqueue;
            
            // init seed points
            for(int label_generator = 0; label_generator < generators.size(); label_generator++)
            {
                Generator& generator = generators[label_generator];
                int point_index = generator.point_index();
                m_vlabels[point_index] = label_generator;
                generator.qem() = m_vqems[point_index];
                add_candidates(pqueue, point_index, label_generator, flag_dist, m_vlabels, generator);
            }

            // partitioning via region growing
            while(!pqueue.empty())
            {
                const CCandidate candidate = pqueue.top();
                pqueue.pop();
                const int point_index = candidate.handle();
                const int label_generator = candidate.index();

                // skip if point already partitioned
                if(m_vlabels.find(point_index) != m_vlabels.end())
                    continue;

                // set label
                m_vlabels[point_index] = label_generator;
                Generator& generator = generators[label_generator];
                generator.add_qem(m_vqems[point_index]);
                add_candidates(pqueue, point_index, label_generator, flag_dist, m_vlabels, Generator);
            }
        }

        /// @brief Add the generators candidates to the priority queue
        /// @param growing_queue 
        /// @param index index of the current point
        /// @param label index of the generator associated to the point
        /// @param flag_dist 
        /// @param m_vlabels 
        /// @param generators 
        void add_candidates(PQueue &pqueue,
        const int index,
        const int label_generator,
        const bool flag_dist,
        const std::map<int, int>& m_vlabels,
        Generator& generator)
        {
            for(const auto neighbor_index : m_graph[index])
            {
                if(m_vlabels.find(neighbor_index) == m_vlabels.end() ) // not already partitioned
                {
                    const double error = compute_growing_error(neighbor_index, generator, flag_dist);
                    pqueue.push(CCandidate(neighbor_index, label_generator, error));
                }
            }
        }

        /// @brief Compute the growing error using the qem cost and weighted by the euclidean distance
        /// @param index index of the current point
        /// @param label index of the generator associated to the point
        /// @param flag flag to use the euclidean distance
        /// @param generators 
        /// @return the cost
        double compute_growing_error(const int neighbor_index,
            Generator& generator,
            const bool flag)
        {
            const double qem_cost = compute_qem_error(m_vqems[neighbor_index], generator.location());
            double total_cost = qem_cost; 

            if(flag)
            {
                Point& neighbor_location = m_point_set.point(neighbor_index);
                total_cost += m_dist_weight * CGAL::squared_distance(generator.location(), neighbor_location);
            }
            return total_cost;
        }

        /// @brief Compute the QEM error from a query point 
        /// @param qem 
        /// @param point 
        /// @return the qem error
        double compute_qem_error(QEM_metric& qem, const Point& point)
        {
            Eigen::VectorXd vec(4);
            vec << point.x(), point.y(), point.z(), 1.0;

            const double error = vec.transpose() * qem.get_4x4_matrix() * vec;
            assert(error >= 0.0);

            return error; 
        }

        /// @brief Update the generators 
        /// @param m_vlabels 
        /// @param generators 
        /// @return a Boolean : true if generators changed and false otherwise 
        bool update_generators(std::map<int, int>& m_vlabels,
            std::vector<Generator>& generators)
        {
            if (m_verbose_level != VERBOSE_LEVEL::LOW)
                std::cout << "Update generators...";

            // records point-indices of current generators
            std::vector<double> min_qem_errors;
            std::vector<int> old_generators;
            for (int label = 0; label < generators.size(); label++)
            {
                Generator& generator = generators[label];
                old_generators.push_back(generator.point_index());
                min_qem_errors.push_back(std::numeric_limits<FT>::max());
            }

            // update generators : optimal qem points for the sum of qem matrices in each cluster
            for (int point_index = 0; point_index < m_point_set.size(); point_index++)
            {
                // skip points not labelled
                if (m_vlabels.find(point_index) == m_vlabels.end())
                    continue;

                // get qem of point's cluster 
                int label = m_vlabels[point_index];
                Generator& generator = generators[label];

                // compute QEM optimal point of generator's cluster, with current point as seed
                Point& seed = m_point_set.point(point_index);
                Point optimal_qem_point = compute_optimal_point(generator.qem(), seed);

                // compute QEM error
                const double qem_error = compute_qem_error(generator.qem(), optimal_qem_point);

                if (qem_error < min_qem_errors[label])
                {
                    generator.point_index() = point_index;
                    generator.location() = optimal_qem_point;
                    min_qem_errors[label] = qem_error;
                }
            }


            // check changes of generators
            for (int i = 0; i < generators.size(); i++)
                if (generators[i].point_index() != old_generators[i])
                    return true;

            if(m_verbose_level > VERBOSE_LEVEL::LOW)
                std::cout << "done" << std::endl;

            // generators have not changed
            return false; 
        }

        // compute errors (total, per cluster, etc)
        void compute_errors()
        {

            // compute errors
            std::vector<double> qem_errors(generators.size(), 0.0);
            for (int i = 0; i < m_point_set.size(); i++)
            {
                if (m_vlabels.find(i) == m_vlabels.end())
                    continue;

                int center_ind = m_vlabels[i];

                double error = compute_minimum_qem_error(m_point_set.point(generators[center_ind]), m_vqems[i]);

                csv_writer->addErrorPoints(i, error);

                if (error > qem_errors[center_ind])
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

        QEM_metric& vqem(const int index)
        {
            assert(index > 0);
            assert(index < m_vqems.size());
            return m_vqems(index);
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