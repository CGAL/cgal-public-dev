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
    Edge(const Node& cell, int c1, int c2, FT v1, FT v2)
    : cell(cell)
    , corner1(c1), corner2(c2)
    , value1(v1), value2(v2)
    , parent(nullptr), child1(nullptr), child2(nullptr) {}

    Edge(const Node& cell, Edge* parent, int c1, int c2, FT v1, FT v2)
    : cell(cell)
    , corner1(c1), corner2(c2)
    , value1(v1), value2(v2)
    , parent(parent), child1(nullptr), child2(nullptr) {}

    ~Edge() {
        delete child1;
        delete child2;
    }

    void divide(const Node& cell1, const Node& cell2, FT midValue) {
        child1 = new Edge(cell1, this, corner1, corner2, value1, midValue);
        child2 = new Edge(cell2, this, corner1, corner2, midValue, value2);
    }

    Edge* findMinimalEdge() {
        if(child1 == nullptr) return this;
        else if(child1->value1 * child1->value2 < 0) return child1->findMinimalEdge();
        else return child2->findMinimalEdge();
    }

    const Node& getCell() const { return cell; }
    Edge* getChild1() const { return child1; }
    Edge* getChild2() const { return child2; }

    std::pair<int,int> corners() const { return std::make_pair(corner1, corner2); }

private:
    const Node& cell;
    int corner1;
    FT value1;
    int corner2;
    FT value2;

    Edge* parent;
    Edge* child1;
    Edge* child2;
};

}

#endif // CGAL_ORTHTREE_EDGE_H