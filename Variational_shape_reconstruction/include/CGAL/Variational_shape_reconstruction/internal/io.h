#ifndef CGAL_VARIATIONAL_SHAPE_RECONSTRUCTION_INTERNAL_IO_H
#define CGAL_VARIATIONAL_SHAPE_RECONSTRUCTION_INTERNAL_IO_H
#include "types.h"

void savePs(std::vector<std::pair<Point, size_t>> m_points,std::map<int, int>   m_vlabels,size_t generator_count,std::string filename)
{
    Pointset pointset;
    std::vector<Vector> colors;
    for(int k = 0 ; k < generator_count;k++)
    {
        colors.push_back(Vector((double) rand() / (RAND_MAX),(double) rand() / (RAND_MAX),(double) rand() / (RAND_MAX)));
    }
    for(int i = 0; i < m_points.size();i++)
    {
        pointset.insert(m_points[i].first,colors[m_vlabels[i]]);
    }
    CGAL::IO::write_XYZ(filename, pointset );
}
    void save_trianglefit_mesh(std::string filename, Polyhedron m_dual_mesh)
    {
        if(m_dual_mesh.empty())
            return;

        std::ofstream mesh_file;
        mesh_file.open(filename);
        CGAL::write_off(mesh_file, m_dual_mesh);
        mesh_file.close();
    }

    void save_trianglefit_soup(std::string filename,std::vector<std::vector<int>> m_facets,PointList m_points)
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

    void save_candidate_edge_ply(std::string filename,std::vector<std::pair<int, int>>m_edges ,PointList m_points)
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
    void save_riemanian(std::string filename,std::vector<std::vector<int>>graph ,PointList m_points)
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