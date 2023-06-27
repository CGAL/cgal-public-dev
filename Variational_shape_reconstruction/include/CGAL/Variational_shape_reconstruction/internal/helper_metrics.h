#include <iostream>
#include <fstream>
#include <vector>
#include <string>

class DataWriter
{
    // iteration 1 iteration 2 .. iteration n
    // pt
    // p2
    // ..
    // pn
    // nb generator
    private:
    std::vector<std::vector<float>> m_data;
    size_t m_points_count;
    public:
    DataWriter(size_t points_count)
    {
        m_points_count = points_count;
        m_data.resize(m_points_count);

    }
    void writeDataToCSV(const std::string& filename) {
        std::ofstream file(filename);

        if (!file) {
            std::cerr << "Error opening the file: " << filename << std::endl;
            return;
        }

        // Write header
        //file << "Name,Age,Email\n";
        // nb of iterions
        size_t nb_iteration = m_data[0].size()-1;
        for(int i = 0 ; i< nb_iteration;i++)
        {
            file << "Iteration "+std::to_string(i)+",";
        }
        file << "Iteration "+std::to_string(nb_iteration)+"\n";

        // Write data
        for (size_t row = 0 ; row < m_points_count;row++) {
            for (size_t line = 0 ; line < nb_iteration ;line++) {
                    file << std::to_string(m_data[row][line])+",";
                }
                file << std::to_string(m_data[row][nb_iteration])+"\n";
        }

        file.close();
    }
    void add (int idx, float error_metric)
    {
        m_data[idx].push_back(error_metric);
    }
};