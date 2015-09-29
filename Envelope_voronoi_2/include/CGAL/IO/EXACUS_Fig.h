#ifndef EXACUS_FIG_H
#define EXACUS_FIG_H

#include <CnX/Conic_point_2.h>
#include <CnX/X_conic_segment_2.h>
#include <CGAL/IO/Fig_stream.h>

#include <vector>

namespace CGAL {

template <class Kernel, class ConicPair_2>
CGAL::Fig_stream<Kernel>& 
  operator<< 
  (
    CGAL::Fig_stream<Kernel>& fig_str, 
    const CnX::Conic_point_2< ConicPair_2 >& point
    )
{
  typedef typename Kernel::Point_2 Point;
  typedef std::pair< double, double > DP;

  if (point.x().is_tending()) 
  {
    return fig_str;
  }    
  if (point.x().number().is_finite()) 
  {
    std::pair< double, double > dp = point.to_double();
    
    // draw only if needed
    if (fig_str.bounding_rect().xmin() <= dp.first &&
        dp.first <= fig_str.bounding_rect().xmax() &&
        fig_str.bounding_rect().ymin() <= dp.second &&
        dp.second <= fig_str.bounding_rect().ymax()) {
      
      fig_str << Point(dp.first, dp.second);
    }
  }
  return fig_str;
}

template <class Kernel, class ConicPoint_2>
CGAL::Fig_stream<Kernel>& 
  operator<< 
  (
    CGAL::Fig_stream<Kernel>& fig_str, 
    const CnX::X_conic_segment_2<ConicPoint_2>& segment
    )
{
  typedef std::pair< double, double > DP;
  
  typedef typename Kernel::Point_2 Point;
  typedef typename Kernel::Segment_2 Segment;

  std::vector< DP > dpoints;
    
  double xmin = fig_str.bounding_rect().xmin();
  double xmax = fig_str.bounding_rect().xmax();
  double ymin = fig_str.bounding_rect().ymin();
  double ymax = fig_str.bounding_rect().ymax();
  double pixel = (xmax - xmin) / fig_str.width();
  segment.discretize_segment(xmin, xmax, ymin, ymax,
                             pixel*10, std::back_inserter(dpoints));
  
  int num_points = static_cast< int >(dpoints.size());
  if (num_points == 0) 
  {
    return fig_str;
  }
    
  if (num_points == 1) 
  {
    // draw single point
    fig_str << segment.source();
  } 
  else
  {
    // draw source and target as point
//    fig_str << segment.source();
//    fig_str << segment.target();
        
    double x0, y0, x1, y1;
    // and connect the point in the list with small line segments
    for (int i = 0;  i + 1 < num_points; i++) {
      x0 = dpoints[i].first;
      y0 = dpoints[i].second;
      x1 = dpoints[i+1].first;
      y1 = dpoints[i+1].second;
      fig_str << Segment(Point(x0,y0), Point(x1,y1));
    }
  }

  return fig_str;
}

} //namespace CGAL

#endif

