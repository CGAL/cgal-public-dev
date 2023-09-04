#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <boost/function_output_iterator.hpp>
#include <boost/property_map/property_map.hpp>
#include <iostream>
#include <fstream>
#include <map>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point;
typedef Kernel::Compare_dihedral_angle_3                    Compare_dihedral_angle_3;
typedef CGAL::Surface_mesh<Point>                           Mesh;
namespace PMP = CGAL::Polygon_mesh_processing;
class DataWriter
{
    // iteration 1 iteration 2 .. iteration n
    // pt
    // p2
    // ..
    // pn
    // nb generator
    private:
    std::vector<std::vector<float>> m_data_error_points;
    std::vector<float> m_data_worst_generators;
    std::vector<float> m_data_mean_generators;
    
    std::map<int,std::vector<float>> m_data_error_generators;
    size_t m_points_count;
    size_t m_generator_count;
    public:
    DataWriter(size_t points_count)
    {
        m_points_count = points_count;
        m_data_error_points.resize(m_points_count);
        //m_data_error_generators.resize(m_points_count);

    }
    /// @brief Compute the mean of the qem error
    /// @param nb_iteration 
    /// @param nb_generator 
    /// @return 
    std::pair<double,double> computeAverage(int nb_iteration, int nb_generator) 
    {
        std::vector<double> values;
        double mean = 0.;
        for (size_t row = 0 ; row < m_generator_count;row++) {
            double v=0.;
            if(m_data_error_generators[row][nb_iteration])
                v = m_data_error_generators[row][nb_iteration];
            if(v<10000000 && v>0)
                values.push_back(v);
        }
        mean = std::accumulate(values.begin(),values.end(),0.);
        mean/=nb_generator;
        // Now calculate the variance
        auto variance_func = [&mean, &nb_generator](double accumulator, const double& val) {
        return accumulator + ((val - mean)*(val - mean) / (nb_generator - 1));
        };
        return {mean,std::accumulate(values.begin(), values.end(), 0.0, variance_func)};
    }
    /// @brief Write the data (qem error for each point) as a csv file
    /// @param filename 
    void writeDataErrorPointsToCSV(const std::string& filename) {
        std::ofstream file(filename);

        if (!file) {
            std::cerr << "Error opening the file: " << filename << std::endl;
            return;
        }

        // Write header
        // nb of iterions
        size_t nb_iteration = m_data_error_points[0].size()-1;
        for(int i = 0 ; i< nb_iteration;i++)
        {
            file << "Iteration "+std::to_string(i)+",";
        }
        file << "Iteration "+std::to_string(nb_iteration)+"\n";

        // Write data
        for (size_t row = 0 ; row < m_points_count;row++) {
            for (size_t line = 0 ; line < nb_iteration ;line++) {
                    file << std::to_string(m_data_error_points[row][line])+",";
                }
                file << std::to_string(m_data_error_points[row][nb_iteration])+"\n";
        }
        //write mean
        file << "nb of generators"+std::to_string(m_generator_count)+"\n mean :,";
        for (size_t line = 0 ; line < nb_iteration ;line++)
        {
            file << std::to_string(computeAverage(line,m_generator_count).first)+",";
            
        }
        file << std::to_string(computeAverage(nb_iteration,m_generator_count).first)+"\n";
        // write variance
        file << "variance:,";
        for (size_t line = 0 ; line < nb_iteration ;line++)
        {
            file << std::to_string(computeAverage(line,m_generator_count).second)+",";
            //std::cout<<"variacne : "<<computeAverage(line,m_generator_count).second<<"\n";
            
        }
        file << std::to_string(computeAverage(nb_iteration,m_generator_count).second)+"\n";
        file.close();
    }
    /// @brief Write the data (qem error for each generator) as a csv file
    /// @param filename 
    void writeDataErrorGeneratorsToCSV(const std::string& filename) {
        std::ofstream file(filename);

        if (!file) {
            std::cerr << "Error opening the file: " << filename << std::endl;
            return;
        }

        // Write header
        // nb of iterions
        size_t nb_iteration =0;
        
        for(auto e : m_data_error_generators)
        {
            nb_iteration = std::max(nb_iteration,e.second.size()-1);
        }
        std::cout<<"nb_iteration "<<nb_iteration<<"\n";
        for(int i = 0 ; i< nb_iteration;i++)
        {
            file << "Iteration "+std::to_string(i)+",";
        }
        file << "Iteration "+std::to_string(nb_iteration)+"\n";

        for(auto e : m_data_error_generators)
        {
            for (auto value : e.second) {
                 file << std::to_string(value)+",";
            }
            file <<"\n";

        }
                //write mean
        file << "nb of generators"+std::to_string(m_generator_count)+"\n mean :,";
        for (size_t line = 0 ; line < nb_iteration ;line++)
        {
            file << std::to_string(computeAverage(line,m_generator_count).first)+",";
            
        }
        file << std::to_string(computeAverage(nb_iteration,m_generator_count).first)+"\n";
        // write variance
        file << "variance:,";
        for (size_t line = 0 ; line < nb_iteration ;line++)
        {
            file << std::to_string(computeAverage(line,m_generator_count).second)+",";
            
        }
        file << std::to_string(computeAverage(nb_iteration,m_generator_count).second)+"\n";
        file.close();
        
    }
    /// @brief print the worst error at each iteration
    void printWorst()
    {
        int id=0;
        for(auto v: m_data_worst_generators)
        {
            std::cout<<"id "<<id<<" error: "<<v<<" mean "<<m_data_mean_generators[id++]<<"\n";
        }
    }
    void addErrorPoints (int idx, float error_metric)
    {
        m_data_error_points[idx].push_back(error_metric);
    }
    void addWorstErrorGenerator (int idx, float error_metric)
    {
        m_data_error_generators[idx].push_back(error_metric);
    }
    void addWorstErrorGenerator (float error_metric)
    {
        m_data_worst_generators.push_back(error_metric);
    }
    void addMeanErrorGenerator (float error_metric)
    {
        m_data_mean_generators.push_back(error_metric);
    }
    void setGenerator(int generator_count)
    {
        m_generator_count= generator_count;
    }
};
