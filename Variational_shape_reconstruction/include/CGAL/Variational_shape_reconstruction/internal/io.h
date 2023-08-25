#ifndef CGAL_VARIATIONAL_SHAPE_RECONSTRUCTION_INTERNAL_IO_H
#define CGAL_VARIATIONAL_SHAPE_RECONSTRUCTION_INTERNAL_IO_H
#include "types.h"

    /// @brief Save the pointset colored by cluster 
    /// @param m_points 
    /// @param m_vlabels 
    /// @param generator_count 
    /// @param filename 
    void savePs(
    const std::vector<std::pair<Point, size_t>>& m_points,
    std::map<int, int>&  m_vlabels,
    const size_t generator_count,
    const std::string filename)
    {
        Pointset pointset;
        std::vector<Vector> colors;
        for(int k = 0 ; k < generator_count;k++)
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
        CGAL::IO::write_XYZ(filename, pointset );
    }
    /// @brief Save the dual mesh 
    /// @param filename 
    /// @param m_dual_mesh 
    void save_trianglefit_mesh(const std::string filename, const Polyhedron& m_dual_mesh)
    {
        if(m_dual_mesh.empty())
            return;

        std::ofstream mesh_file;
        mesh_file.open(filename);
        CGAL::write_off(mesh_file, m_dual_mesh);
        mesh_file.close();
    }
    /// @brief Save the soup of triangle of m_facets
    /// @param filename 
    /// @param m_facets list of the indices of each facet
    /// @param m_points 
    void save_trianglefit_soup(const std::string filename,
    const std::vector<std::vector<int>>& m_facets,
    const PointList& m_points)
    {
        if(m_facets.size() == 0)
            return;

        std::ofstream soup_file;
        soup_file.open(filename);

        soup_file << "OFF"  << std::endl;
        soup_file << m_points.size() << " " << m_facets.size() << " 0" << std::endl;

        for(int i = 0; i < m_points.size(); i++)
            soup_file << m_points[i].x() << " " << m_points[i].y() << " " << m_points[i].z() << std::endl;

        for(int i = 0; i < m_facets.size(); i++)
            soup_file << "3 " << m_facets[i][0] << " " << m_facets[i][1] << " " << m_facets[i][2] << std::endl;

        soup_file.close();
    }
    /// @brief save the candidate edge as ply
    /// @param filename 
    /// @param m_edges 
    /// @param m_points 
    void save_candidate_edge_ply(const std::string filename,
    const std::vector<std::pair<int, int>>& m_edges,
    const PointList& m_points)
    {
        std::ofstream edge_file;
        edge_file.open(filename);

        edge_file << "ply\n"
                  << "format ascii 1.0\n"
                  << "element vertex " << m_points.size() << "\n"
                  << "property float x\n"
                  << "property float y\n"
                  << "property float z\n"
                  << "element edge " << m_edges.size() << "\n"
                  << "property int32 vertex1\n"
                  << "property int32 vertex2\n"
                  << "end_header\n";

        for(int i = 0; i < m_points.size(); i++)
            edge_file << m_points[i].x() << " " << m_points[i].y() << " " << m_points[i].z() << std::endl;
        
        for(int i = 0; i < m_edges.size(); i++)
            edge_file << m_edges[i].first << " " << m_edges[i].second << std::endl;

        edge_file.close();
    }
    /// @brief Save the rimemanian graph 
    /// IE the graph of the neighbors as edge between the point and its neighbors
    /// @param filename 
    /// @param graph 
    /// @param m_points 
    void save_riemanian(const std::string filename,
    copnst std::vector<std::vector<int>>& graph,
    const PointList& m_points)
    {
        std::ofstream edge_file;
        edge_file.open(filename);

        std::size_t sum = 0;
        for (auto &&i : graph) {
            sum += i.size();
        }

        edge_file << "ply\n"
                  << "format ascii 1.0\n"
                  << "element vertex " << m_points.size() << "\n"
                  << "property float x\n"
                  << "property float y\n"
                  << "property float z\n"
                  << "element face " << sum << "\n"
                  << "property list uchar int vertex_indices\n"
                  << "end_header\n";

        for(int i = 0; i < m_points.size(); i++)
            edge_file << m_points[i].x() << " " << m_points[i].y() << " " << m_points[i].z() << std::endl;
        
        for(int i = 0; i < graph.size(); i++)
        {
            for(int j = 0; j < graph[i].size(); j++)
            {
                edge_file << "2 "<<i<<" "<<j<<"\n";
            }
        }

        edge_file.close();
    }

#endif // CGAL_VARIATIONAL_SHAPE_RECONSTRUCTION_INTERNAL_IO_H