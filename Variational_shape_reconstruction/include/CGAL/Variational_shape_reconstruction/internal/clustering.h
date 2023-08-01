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
class Clustering
{
    public:
    Clustering()
    {
    }
    Clustering(const Pointset& pointset, size_t num_knn, double euclidean_distance_weight)
    {
        pointset_ = pointset;
        //m_num_knn = num_knn;
        m_dist_weight = euclidean_distance_weight;
        //csv_writer =std::make_shared<DataWriter>(pointset.size());
    }
    void initialize_qem_map(const KNNTree& m_tree)
    {
        int num_nb = std::max(6, m_num_knn);
        for(Pointset::const_iterator it = pointset_.begin(); it != pointset_.end(); ++ it)
        {
            auto point = pointset_.point(*it);
            K_neighbor_search search(m_tree, point, num_nb);
            KNNDistance tr_dist;

            double avg_dist = 0.;
            for(typename K_neighbor_search::iterator it = search.begin(); it != search.end(); it++)
                avg_dist += tr_dist.inverse_of_transformed_distance(it->second);

            avg_dist = avg_dist / (double)num_nb;
            
            QEM_metric pqem = compute_qem_for_point(point, pointset_.normal(*it), avg_dist * avg_dist);
            m_pqems.push_back(pqem);
        }
    }
    void compute_connected()
    {
        size_t point_count = m_graph.size(); 
        m_visited.resize(point_count);
        m_component.resize(point_count);
        std::fill(m_visited.begin(),m_visited.end(),false);
        
        for(int i = 0 ; i < point_count;i++)
        {
            if(!m_visited[i])
            {
                visit(i,component_count);
                component_count++;
            }
        }
    }
    void visit(int current_idx, int component_idx)
    {
        std::vector<int> queue;
        queue.push_back(current_idx);
        
        while(!queue.empty())
        {
            int idx = queue.back();
            queue.pop_back();
            m_component[idx]=component_idx;
            m_visited[idx]=true;
            for(auto neighbors_idx : m_graph[idx])
            {
                if(!m_visited[neighbors_idx])
                {
                    queue.push_back(neighbors_idx);
                }
            }
        }
    }
    QEM_metric compute_qem_for_point(const Point& query,const Vector& normal,const double &area)
    {
        QEM_metric qem;
        qem.init_qem_metrics_face(area, query, normal);
        return qem;
    }
    void initialize_vertex_qem(const KNNTree& m_tree,std::vector<int>& m_generators)
    {
        m_vqems.clear();
        for(Pointset::const_iterator it = pointset_.begin(); it != pointset_.end(); ++ it)
        {
            K_neighbor_search search(m_tree, pointset_.point(*it), m_num_knn);
            QEM_metric vqem;
            std::vector<int> neighbors;
            for(typename K_neighbor_search::iterator it = search.begin(); it != search.end(); it++)
            {
                auto neighbor_idx = (it->first).second;
                vqem = vqem + m_pqems[neighbor_idx];
                neighbors.push_back(neighbor_idx);
            }
            m_graph.push_back(neighbors);

            m_vqems.push_back(vqem);
        }
        compute_connected();
        std::ofstream edge_file;
        edge_file.open("filename.ply");

        std::size_t sum = 0;
        for (auto &&i : m_graph) {
            sum += i.size();
        }

        edge_file << "ply\n"
                  << "format ascii 1.0\n"
                  << "element vertex " << pointset_.size() << "\n"
                  << "property float x\n"
                  << "property float y\n"
                  << "property float z\n"
                  << "element face " << sum << "\n"
                  << "property list uchar int vertex_indices\n"
                  << "end_header\n";

        for(Pointset::const_iterator it = pointset_.begin(); it != pointset_.end(); ++ it)
        {
            auto point = pointset_.point(*it);
            edge_file << point.x() << " " << point.y() << " " << point.z() << std::endl;
        }
        
        for(int i = 0; i < m_graph.size(); i++)
        {
            for(int j = 0; j < m_graph[i].size(); j++)
            {
                edge_file << "2 "<<i<<" "<<m_graph[i][j]<<"\n";
            }
        }

        edge_file.close();
        
        std::cout<<"connected component count "<< component_count<<std::endl;
        for(int i = 0 ; i < component_count;i++)
        {
            bool found = false;
            for(int j = 0 ; j < m_generators.size();j++)
            {
                if(i == m_component[m_generators[j]])
                {
                    found = true;
                }
            }
            if(!found)
            {
                // Add a new generator in the connected component without generator
                auto it = std::find(m_component.begin(), m_component.end(), i);
                int index = it - m_component.begin();
                m_generators.push_back(index);
                std::cout<<"added generator > "<<index<<" connected component "<< i<<std::endl;
            }
        }
    }
    void region_growing(std::map<int, int>& m_vlabels,std::vector<QEM_metric>& m_generators_qem,const std::vector<int>& m_generators,bool flag_dist)
    {
        PQueue growing_queue;
        
        // init seed triangles
        for(int label = 0; label < m_generators.size(); label++)
        {
            int index = m_generators[label];
            m_vlabels[index] =  label;
            m_generators_qem.push_back(m_vqems[index]);
            add_candidates(growing_queue,index,label, flag_dist,m_vlabels,m_generators);
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
            m_generators_qem[label] = m_generators_qem[label] + m_vqems[index]; 
            add_candidates(growing_queue,index,label, flag_dist,m_vlabels,m_generators);
        }


    }
    void add_candidates(PQueue &growing_queue,const int index,const int label,const bool flag_dist,
                        const std::map<int, int>& m_vlabels,const std::vector<int>& m_generators)
    {
        
        for(const auto nb_index : m_graph[index])
        {
            if(m_vlabels.find(nb_index) == m_vlabels.end() )
            {
                const double loss = compute_collapse_loss(nb_index, label, flag_dist,m_generators);
                growing_queue.push(CCandidate(nb_index, label, loss));
            }
        }
    }
    double compute_collapse_loss(const int index,const int label,const bool flag,const std::vector<int>& m_generators) //eq 6
    {
        const double m_qem_weight=1.;
        
        
        const double qem_cost = compute_minimum_qem_error(pointset_.point(m_generators[label]), m_vqems[index]);

         double cost = m_qem_weight * qem_cost; 

        if(flag)
        {
            double dist_cost = m_num_knn * CGAL::squared_distance(pointset_.point(m_generators[label]), pointset_.point(index));
            cost = cost + m_dist_weight * dist_cost;
        }
        return cost;
    }
    double compute_minimum_qem_error(Point center_point, QEM_metric& query_qem)
    {
        Eigen::VectorXd center_vec(4);
        center_vec << center_point.x(), center_point.y(), center_point.z(), 1.;

        double cost = center_vec.transpose() * query_qem.get_4x4_matrix() * center_vec;

        return std::abs(cost);
    }
    bool update_poles(std::map<int, int>& m_vlabels,std::vector<QEM_metric>& m_generators_qem,std::vector<int>& m_generators)
    {
        std::cout << "Updating poles..." << std::endl;

        std::vector<Point>    optimal_points;
        std::vector<double> dists;
        std::vector<int> old_poles;

        for(int i = 0; i < m_generators.size(); i++)
        {
            Point center;
            // is it for multicover ?
            /*if(m_generators_count.size() == m_generators.size() && m_generators_count[i] == 1)
            center = m_points[m_generators[i]].first;
            else*/
            center = compute_optimal_point(m_generators_qem[i], pointset_.point(m_generators[i]));


            //center = pointset_.point(m_generators[i]);
            optimal_points.push_back(center);
            dists.push_back(1e20);
            old_poles.push_back(m_generators[i]);
        }
        
        for(int i = 0; i < pointset_.size(); i++) 
        {
            if(m_vlabels.find(i) == m_vlabels.end())
                continue;

            int label = m_vlabels[i];
            double dist = CGAL::squared_distance(optimal_points[label], pointset_.point(i));

            if(dist < dists[label])
            {
                m_generators[label] = i;
                dists[label] = dist;
            }
        }
        /// analysis
        /*for(int i = 0; i < pointset_.size(); i++) 
        {
            if(m_vlabels.find(i) == m_vlabels.end())
                continue;

            int center_ind = m_vlabels[i];
            double error = compute_minimum_qem_error(pointset_.point(m_generators[center_ind]), m_vqems[i]); 
            //csv_writer->addErrorPoints(i,error);

        }*/
        std::vector<double> qem_errors(m_generators.size(), 0.);
        for(int i = 0; i < pointset_.size(); i++) 
        {
            if(m_vlabels.find(i) == m_vlabels.end())
                continue;

            int center_ind = m_vlabels[i];
            double error = compute_minimum_qem_error(pointset_.point(m_generators[center_ind]), m_vqems[i]); 
            //csv_writer->addErrorPoints(i,error);

            if(error > qem_errors[center_ind])
            {
                qem_errors[center_ind] = error;
            }
        }
        for(int i = 0 ; i < qem_errors.size();i++)
        {
            //csv_writer->addWorstErrorGenerator(i,qem_errors[i]);
        }
        std::cout<<"Generators: "<<m_generators.size()<<"\n";
        //csv_writer->setGenerator(m_generators.size());
        // chech change
        for(int i = 0; i < m_generators.size(); i++)
        {
            if(m_generators[i] != old_poles[i])
                return true;
        }

        //std::cout << "Region growing converges!" << std::endl;
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
    size_t guided_split_clusters(/*const Pointset& pointset,const std::vector<int> generators,const std::vector<QEM_metric> vqems,*/
                            std::map<int, int>& m_vlabels,std::vector<int>& m_generators,
                             double m_diag,double m_spacing,
                             double split_ratio, size_t iteration) // batch splitting
{


    double split_thresh = m_diag * split_ratio;
    split_thresh = split_thresh * split_thresh; // square distance
    split_thresh = split_thresh * m_num_knn * m_spacing * m_spacing;
    //std::cout << "  split threshold: " << split_thresh << std::endl;

    // compute error for each center
    std::vector<double> qem_errors(m_generators.size(), 0.);
    std::vector<int> vert_indices(m_generators.size(), -1);

    std::map<int,double> generator_worst_error;
    for(int i = 0; i < pointset_.size(); i++) 
    {
        if(m_vlabels.find(i) == m_vlabels.end())
            continue;

        int center_ind = m_vlabels[i];
        double error = compute_minimum_qem_error(pointset_.point(m_generators[center_ind]), m_vqems[i]); 
        //csv_writer->addErrorPoints(i,error);
        /*if(generator_worst_error.count(center_ind) > 0)
            generator_worst_error[center_ind]= std::min(generator_worst_error[center_ind],error);
        else
            generator_worst_error[center_ind] = error;
            */


        if(error > qem_errors[center_ind])
        {
            qem_errors[center_ind] = error;
            vert_indices[center_ind] = i;
        }
    }
    
    for(int i = 0 ; i < qem_errors.size();i++)
    {
        //csv_writer->addWorstErrorGenerator(i,qem_errors[i]);
    }
    std::cout<<"Generators: "<<m_generators.size()<<"\n";
    //csv_writer->setGenerator(m_generators.size());

    // split centers exceeding max error
    std::vector<int> new_poles;

    for(int i = 0; i < m_generators.size(); i++)
    {

        if(qem_errors[i] > split_thresh && vert_indices[i] != m_generators[i])
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
                CGAL::squared_distance(pointset_.point(new_poles[i]), pointset_.point(new_poles[j])) < dist_thresh)
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
        m_generators.push_back(new_poles[i]);

    }

    std::cout << "Finally found " << new_poles.size() << " new poles!" << std::endl;


    return new_poles.size();
}
void write_csv()
{
    //csv_writer->writeDataErrorGeneratorsToCSV("error_generators.csv");
    //csv_writer->writeDataErrorPointsToCSV("error_points.csv");
}

        private:
            Pointset pointset_;
            int m_num_knn = 12;
            double m_dist_weight=0.1;
                    //qem
            std::vector<QEM_metric> m_pqems;
            std::vector<QEM_metric> m_vqems;
            std::vector<std::vector<int>> m_graph;
            std::vector<bool> m_visited;
            std::vector<int> m_component;
            int component_count=0;
          // csv

            //std::shared_ptr<DataWriter> csv_writer;
            
            


};