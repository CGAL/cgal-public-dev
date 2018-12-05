// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYLINE_TRACING_MG_VALIDITY_CHECK_H
#define CGAL_POLYLINE_TRACING_MG_VALIDITY_CHECK_H

#include <CGAL/assertions.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/optional.hpp>

#include <iostream>

namespace CGAL {

namespace Polyline_tracing {

namespace internal {

// mega-brute-force spaghetti-code validity check

template<typename MotorcycleGraph>
bool is_valid(const MotorcycleGraph& mg)
{
  typedef MotorcycleGraph                                                   Motorcycle_graph;

  typedef typename Motorcycle_graph::Geom_traits                            Geom_traits;
  typedef typename Geom_traits::FT                                          FT;
  typedef typename Geom_traits::Point_2                                     Point_2;
  typedef typename Geom_traits::Segment_2                                   Segment_2;

  typedef typename Geom_traits::Triangle_mesh                               Triangle_mesh;
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;

  typedef typename Motorcycle_graph::Node_ptr                               Node_ptr;
  typedef typename Motorcycle_graph::Motorcycle                             Motorcycle;
  typedef typename Motorcycle_graph::Track                                  Track;
  typedef typename Motorcycle_graph::Track_segment                          Track_segment;

  typedef typename Motorcycle_graph::MCC_cit                                MCC_cit;

  namespace PMP = CGAL::Polygon_mesh_processing;

  std::cout << "Checking trace validity..." << std::endl;

  MCC_cit mc_it = mg.motorcycles().begin(), mc_end = mg.motorcycles().end();
  for(; mc_it!=mc_end; ++mc_it)
  {
    const Motorcycle& mc = mg.motorcycle(mc_it);
    const Track& mc_track = mc.track();
    std::size_t mc_track_size = mc_track.size();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "Track of motorcycle: " << mc.id() << std::endl << mc_track << std::endl;
#endif

    assert(mc_track.is_valid());

    if(mc_track_size == 0)
      continue;

    // If the motorcycle did not reach its destination, it must have crashed on the track of another motorcycle
    if(mc.current_position() != mc.destination())
      CGAL_assertion(mc.current_position()->earliest_visiting_time() <= mc.current_time());

    typename Track::const_iterator mc_tcit = mc_track.begin(), end = mc_track.end();
    for(; mc_tcit!=end; ++mc_tcit)
    {
      const Track_segment& mc_ts = *mc_tcit;

      const Node_ptr mc_track_source = mc_ts.source();
      const Node_ptr mc_track_target = mc_ts.target();
      const FT time_at_mc_track_target = mc_ts.time_at_target();
      const face_descriptor fd = mc_track_source->face();

      // Ignore degenerate tracks, except if the whole track is a single degenerate segment
      const bool is_mc_degenerate = (mc_track_source == mc_track_target);
      if(is_mc_degenerate && mc_track_size > 1)
        continue;

      // Here comes the mega brute force: check each track segment against ALL the other track segments
      const Point_2 ts = mg.geom_traits().construct_point_2_object()(mc_track_source->barycentric_coordinate(0),
                                                                     mc_track_source->barycentric_coordinate(1));
      const Point_2 tt = mg.geom_traits().construct_point_2_object()(mc_track_target->barycentric_coordinate(0),
                                                                     mc_track_target->barycentric_coordinate(1));
      Segment_2 s = mg.geom_traits().construct_segment_2_object()(ts, tt);

      MCC_cit fmc_it = mg.motorcycles().begin(), fmc_end = mg.motorcycles().end();
      for(; fmc_it!=fmc_end; ++fmc_it)
      {
        const Motorcycle& fmc = mg.motorcycle(fmc_it);
        const Track& fmc_track = fmc.track();
        std::size_t fmc_track_size = fmc_track.size();

        if(fmc_track_size == 0)
          continue;

        typename Track::const_iterator fmc_tcit = fmc_track.begin(), fend = fmc_track.end();
        for(; fmc_tcit!=fend; ++fmc_tcit)
        {
          const Track_segment& fmc_ts = *fmc_tcit;

          const Node_ptr fmc_track_source = fmc_ts.source();
          const FT time_at_fmc_track_source = fmc_ts.time_at_source();
          const Node_ptr fmc_track_target = fmc_ts.target();
          const FT time_at_fmc_track_target = fmc_ts.time_at_target();
          const face_descriptor ffd = fmc_track_source->face();

          if(mc_ts == fmc_ts)
            continue;

          // Ignore degenerate segments, except if the whole track is a single degenerate segment
          bool is_fmc_degenerate = (fmc_track_source == fmc_track_target);
          if(is_fmc_degenerate && fmc_track_size > 1)
            continue;

          Segment_2 fs;
          if(fd == ffd)
          {
            const Point_2 fts = mg.geom_traits().construct_point_2_object()(fmc_track_source->barycentric_coordinate(0),
                                                                            fmc_track_source->barycentric_coordinate(1));
            const Point_2 ftt = mg.geom_traits().construct_point_2_object()(fmc_track_target->barycentric_coordinate(0),
                                                                            fmc_track_target->barycentric_coordinate(1));
            fs = mg.geom_traits().construct_segment_2_object()(fts, ftt);
          }
          else // on different faces
          {
            boost::optional<halfedge_descriptor> ohd = PMP::common_halfedge(fd, ffd, mg.mesh());
            if(ohd != boost::none) // the two faces are adjacent
            {
              halfedge_descriptor chd = *ohd;
              halfedge_descriptor ochd = opposite(chd, mg.mesh());

              bool is_source_on_common_edge = PMP::is_on_halfedge(mc_track_source->location(), chd, mg.mesh());
              bool is_target_on_common_edge = PMP::is_on_halfedge(mc_track_target->location(), chd, mg.mesh());

              if(!is_source_on_common_edge && !is_target_on_common_edge)
                continue;

              bool is_fsource_on_common_edge = PMP::is_on_halfedge(fmc_track_source->location(), ochd, mg.mesh());
              bool is_ftarget_on_common_edge = PMP::is_on_halfedge(fmc_track_target->location(), ochd, mg.mesh());

              std::cout << "foreign track segment on foreign but adjacent face " << ffd << std::endl;
              std::cout << "on common halfedge: " << is_fsource_on_common_edge << " " << is_ftarget_on_common_edge << std::endl;

              if(is_fsource_on_common_edge)
              {
                if(is_ftarget_on_common_edge)
                {
                  // foreign track completely on the edge
                  const Node_ptr fsource_in_fd = fmc_track_source->sibling(fd);
                  const Node_ptr ftarget_in_fd = fmc_track_target->sibling(fd);

                  const Point_2 fts = mg.geom_traits().construct_point_2_object()(fsource_in_fd->barycentric_coordinate(0),
                                                                                  fsource_in_fd->barycentric_coordinate(1));
                  const Point_2 ftt = mg.geom_traits().construct_point_2_object()(ftarget_in_fd->barycentric_coordinate(0),
                                                                                  ftarget_in_fd->barycentric_coordinate(1));
                  fs = mg.geom_traits().construct_segment_2_object()(fts, ftt);
                }
                else // !is_ftarget_on_common_edge
                {
                  // only the source is on the common halfedge, create a degenerate track segment
                  const Node_ptr fsource_in_fd = fmc_track_source->sibling(fd);

                  const Point_2 fts = mg.geom_traits().construct_point_2_object()(fsource_in_fd->barycentric_coordinate(0),
                                                                                  fsource_in_fd->barycentric_coordinate(1));
                  fs = mg.geom_traits().construct_segment_2_object()(fts, fts);
                  is_fmc_degenerate = true;
                }
              }
              else // !is_fsource_on_common_edge
              {
                // only the target is on the common halfedge, create a degenerate track segment
                if(is_ftarget_on_common_edge)
                {
                  const Node_ptr ftarget_in_fd = fmc_track_target->sibling(fd);

                  const Point_2 ftt = mg.geom_traits().construct_point_2_object()(ftarget_in_fd->barycentric_coordinate(0),
                                                                                  ftarget_in_fd->barycentric_coordinate(1));
                  fs = mg.geom_traits().construct_segment_2_object()(ftt, ftt);
                  is_fmc_degenerate = true;
                }
                else // !is_ftarget_on_common_edge
                {
                  // foreign track segment does not touch the border of the current track's face
                  continue;
                }
              }

            }
            else // faces are not adjacent, no chance of collision
            {
              continue;
            }
          }

          // Actually check for collision between the track segments
          if(mg.geom_traits().do_intersect_2_object()(s, fs))
          {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
            std::cout << "Intersection ¤~~~~~~~~~~~~~~~~~¤" << std::endl;
            std::cout << "motorcycle #" << mc.id() << " (track size: " << mc_track.size();
            std::cout << ") with motorcycle #" << fmc.id() << " (track size: " << fmc_track.size() << ")" << std::endl;
            std::cout << "DECITs:" << std::endl
                      << "-------" << std::endl << *mc_track_source << std::endl
                      << "-------" << std::endl << *mc_track_target << std::endl
                      << "-------" << std::endl << *fmc_track_source << std::endl
                      << "-------" << std::endl << *fmc_track_target << std::endl;
#endif

            // An intersection must be at an extremity
            if(is_mc_degenerate)
            {
              if(is_fmc_degenerate) // both track segments are generate
              {
                CGAL_assertion(s.source() == fs.source());
              }
              else // only mc's track segment is denegerate
              {
                CGAL_assertion((fs.source() == s.source() && fs.target() != s.source()) ||
                               (fs.source() != s.source() && fs.target() == s.source()));
              }
            }
            else if(is_fmc_degenerate) // only the foreign track segment is denegerate
            {
              CGAL_assertion((s.source() == fs.source() && s.target() != fs.source()) ||
                             (s.source() != fs.source() && s.target() == fs.source()));
            }
            else // no degenerate track segment
            {
              CGAL_assertion((s.source() == fs.source() && s.source() != fs.target() &&
                              s.target() != fs.source() && s.target() != fs.target()) ||
                             (s.source() != fs.source() && s.source() == fs.target() &&
                              s.target() != fs.source() && s.target() != fs.target()) ||
                             (s.source() != fs.source() && s.source() != fs.target() &&
                              s.target() == fs.source() && s.target() != fs.target()) ||
                             (s.source() != fs.source() && s.source() != fs.target() &&
                              s.target() != fs.source() && s.target() == fs.target()));
            }

            // An intersection at the target of the track must crash the motorcycle
            // if the motorcycle reaches this collision point at a later time
            // than another motorcycle.
            // In other words, if there is an intersection at 'mc_track_target'
            // AND the visiting time is lower for the foreign motorcycle, then
            // 'mc_track_target' must be the last track entry for 'mc'.

            bool is_fmc_track_segment_next_after_current = (mc.id() == fmc.id() &&
                                                            s.target() == fs.source());

            if(!is_fmc_track_segment_next_after_current) // otherwise false positives
            {
              // only '<' for the foreign source because equal times are acceptable
              if((mc_track_target->point() == fmc_track_source->point() &&
                  time_at_fmc_track_source < time_at_mc_track_target) ||
                 (mc_track_target->point() == fmc_track_target->point() &&
                   (time_at_fmc_track_target < time_at_mc_track_target ||
                     (time_at_fmc_track_target == time_at_mc_track_target && !is_fmc_degenerate))))
              {
                // must be the last item of the track
                if(mc_tcit != --(mc_track.end()))
                {
                  std::cerr << "Motorcycle: " << std::endl << mc << std::endl;
                  std::cerr << "should have been stopped at: " << std::endl << *mc_track_target << std::endl;
                  std::cerr << "by foreign motorcycle : " << std::endl << fmc << std::endl;
                  if(mc_track_target->point() == fmc_track_source->point())
                    std::cerr << "times: " << time_at_mc_track_target << " vs " << time_at_fmc_track_source << std::endl;
                  else
                    std::cerr << "times: " << time_at_mc_track_target << " vs " << time_at_fmc_track_target << std::endl;

                  CGAL_assertion(false);
                }
              }
            }
          }
        }
      }
    }
  }

  CGAL::Polyline_tracing::internal::is_valid_hds(mg.graph());

  return true;
}


} // namespace internal

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MG_VALIDITY_CHECK_H
