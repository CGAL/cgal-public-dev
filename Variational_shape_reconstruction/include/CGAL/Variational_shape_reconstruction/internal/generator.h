#ifndef GENERATOR_H
#define GENERATOR_H

#include "types.h"
#include "qem.h"

namespace qem
{
    template <class Kernel>
    class CGenerator
    {
    private:
        int m_index; //  index in input point set
        Point m_location; // optimal location in space
        QEM_metric m_qem; // sum of qem of the whole cluster
    public:
        int& index() { return m_index; }
        Point& location() { return m_location; }
        QEM_metric& qem() { return m_qem; }

        // default constructor
        CGenerator()
        {
            m_index = 0;
            m_location = CGAL::ORIGIN;
        }

        // constructor with null qem
        CGenerator(const int index,
            const Point& location)
        {
            m_index = index;
            m_location = location;
        }

        // copy constructor
        CGenerator(CGenerator& g)
            : m_index(g.index()), 
            m_location(g.location()), 
            m_qem(g.qem())
        {
        }
    };
}

#endif /* GENERATOR_H */
