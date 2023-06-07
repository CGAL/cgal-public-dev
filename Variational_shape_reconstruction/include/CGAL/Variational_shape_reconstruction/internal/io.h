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

#endif // CGAL_VARIATIONAL_SHAPE_RECONSTRUCTION_INTERNAL_IO_H