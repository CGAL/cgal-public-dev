#ifndef CGAL_ORTHTREE_EDGE_H
#define CGAL_ORTHTREE_EDGE_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Orthtree.h>

namespace CGAL {

/*!

  \brief represents an edge of the tree.
 */
template<typename Traits, typename PointRange, typename PointMap>
class Orthtree<Traits, PointRange, PointMap>::Edge
{
    typedef typename Traits::FT FT;
    typedef CGAL::Vector_3<Traits> Vector;

public:
    Edge(const std::array<std::uint32_t, 6>& coords, int depth, FT v1, FT v2)
    : Global_coordinates(coords), Depth(depth)
    , value1(v1), value2(v2)
    , parent(nullptr), child1(nullptr), child2(nullptr) {}

    Edge(const std::array<std::uint32_t, 6>& coords, int depth, Edge* parent, FT v1, FT v2)
    : Global_coordinates(coords), Depth(depth)
    , value1(v1), value2(v2)
    , parent(parent), child1(nullptr), child2(nullptr) {}

    ~Edge() {
        delete child1;
        delete child2;
    }

    void divide(FT midValue) {
        std::array<std::uint32_t, 3> midCoords { (Global_coordinates[0] + Global_coordinates[3]),
                    (Global_coordinates[1] + Global_coordinates[4]),
                    (Global_coordinates[2] + Global_coordinates[5]) };
        child1 = new Edge(
            {2*Global_coordinates[0], 2*Global_coordinates[1], 2*Global_coordinates[2], midCoords[0], midCoords[1], midCoords[2]}
            , Depth + 1, this, value1, midValue);
        child2 = new Edge(
            {midCoords[0], midCoords[1], midCoords[2], 2*Global_coordinates[3], 2*Global_coordinates[4], 2*Global_coordinates[5]}
            , Depth + 1, this, midValue, value2);
    }

    Edge* findMinimalEdge() {
        if(child1 == nullptr) return this;
        else if(child1->value1 * child1->value2 <= 0) return child1->findMinimalEdge();
        else return child2->findMinimalEdge();
    }

    Edge* getChild1() const { return child1; }
    Edge* getChild2() const { return child2; }

    // the two corner points of the given node that are on this edge
    std::pair<int,int> corners(Node n) const {
        int corner1 = (Global_coordinates[0] > n.global_coordinates()[0])*4 
            + (Global_coordinates[1] > n.global_coordinates()[1])*2
            + (Global_coordinates[2] > n.global_coordinates()[2])*1;
        int corner2 = (Global_coordinates[3] > n.global_coordinates()[0])*4 
            + (Global_coordinates[4] > n.global_coordinates()[1])*2
            + (Global_coordinates[5] > n.global_coordinates()[2])*1;
        return std::make_pair(corner1, corner2);
    }
    std::array<std::uint32_t, 6> global_coordinates() const { return Global_coordinates; }
    int depth() const { return Depth; }

    friend std::ostream& operator<<(std::ostream& out, const Edge& edge) {
        auto c = edge.global_coordinates();
        out << "{ ( " << c[0] << " " << c[1] << " " << c[2] << " ) ( " << c[3] << " " << c[4] << " " << c[5] << " ) " << edge.depth() << " }" << std::endl;
        return out;
    }

private:
    std::array<std::uint32_t, 6> Global_coordinates;
    int Depth;
    FT value1;
    FT value2;

    Edge* parent;
    Edge* child1;
    Edge* child2;
};

}

#endif // CGAL_ORTHTREE_EDGE_H