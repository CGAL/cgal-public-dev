// Offset_statistics_2.h

#ifndef Offset_statistics_2_h
#define Offset_statistics_2_h

#include <CGAL/Offset_types_2.h>
#include <CGAL/Timer.h>
#include <iostream>

namespace CGAL {


/*! \classf
 * A class implementing statistics govering for various polygons
 * and polygon offsets computations. 
 */

template <typename CPolygon_2>
class Offset_statistics_2 : public Offset_types_2<CPolygon_2>
{
public:

  typedef Offset_types_2<CPolygon_2> Types;

  Offset_statistics_2(void) : m_verbose(true)
  {
  }

  ~Offset_statistics_2(void)
  {
  }

  struct Polygon_statistics
  {
      // computation time
      double time;
      // number of outer boundary vertices
      unsigned int size;

      Polygon_statistics() : time(0.0), size(0) {}
  };

  struct Offset_statistics
  {
      // computation time
      double time;
      // number of outer boundary vertices
      unsigned int size;
      // number of holes
      unsigned int holes;
      // number of arcs
      unsigned int circles;

      Offset_statistics() : time(0.0), size(0), holes(0), circles(0) {}
  };

  virtual void verbose(const bool& be_verbose) { m_verbose = be_verbose; }
  const bool& verbose() const { return m_verbose; }

protected:
  
  bool m_verbose; // true if should output statistics after performed computations

private:

};

} // namespace CGAL


#endif // Offset_statistics_2_h
