#ifndef CGAL_LTS_INTERNAL_H
#define CGAL_LTS_INTERNAL_H

#include <utility>

#include <boost/array.hpp>

namespace CGAL { namespace LTS {

template <typename OutputIterator, typename Transversal, typename Segment>
OutputIterator insert_transversal(
  OutputIterator out, const Transversal& output_segment,
  const Segment* s1, const Segment* s2, 
  const Segment* s3, const Segment* s4,
  boost::true_type) 
{
  boost::array<const Segment*, 4> arr = {{s1, s2, s3, s4}};
  return *out++ = std::make_pair(output_segment, arr);
}

template <typename OutputIterator, typename Transversal, typename Segment>
OutputIterator insert_transversal(
  OutputIterator out, const Transversal& output_segment,
  const Segment* s1, const Segment* s2, 
  const Segment* s3, const Segment* s4,
  boost::false_type) 
{
  return *out++ = output_segment;
}


} // LTS
} // CGAL


#endif /* CGAL_LTS_INTERNAL_H */
