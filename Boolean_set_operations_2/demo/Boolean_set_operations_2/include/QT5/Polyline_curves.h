// Copyright (c) 2012  Tel-Aviv University (Israel).
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
// Author(s) : Apurva Bhatt <response2apurva@gmail.com>
//             Ronnie Gandhi <ronniegandhi19999@gmail.com>
//             Efi Fogel <efifogel@gmain.com>

#pragma once


#include <CGAL/Qt/Converter.h>

#include "Typedefs.h"
#include "QT5/Piecewise_set_graphics_item.h"
#include "QT5/Boundary_pieces_graphics_item.h"

namespace CGAL {
    namespace Qt {

        struct Polyline_X_monotone_bbox {
            template <typename X_monotone_polyline_segment_2>
            CGAL::Bbox_2 operator()(X_monotone_polyline_segment_2 const& aC) const
            { return aC.bbox(); }
        };

        struct Polyline_bbox {
            template<class Polyline_segment_2>
            CGAL::Bbox_2 operator()(Polyline_segment_2 const& aC) const
            {
              /*  double x_min = to_double(aC.source().x());
                double x_max = to_double(aC.target().x());
                double y_min = to_double(aC.source().y());
                double y_max = to_double(aC.target().y());
                if(x_min > x_max) {
                    std::swap(x_min, x_max);
                    std::swap(y_min, y_max);
                }

                if(y_min > y_max) std::swap(y_min, y_max);

                return Bbox_2(x_min, y_min, x_max, y_max);*/
                auto n = aC.number_of_subcurves();
                Bbox_2 bbox;
                for (std::size_t i = 0; i < n; ++i)
                    bbox = (i > 0) ? (bbox + aC[i].bbox()) : aC[i].bbox();
                return bbox;
            }
        };

        struct Draw_polyline_X_monotone_curve {
            template <typename X_monotone_polyline_segment_2, typename Path>
            void operator()(X_monotone_polyline_segment_2 const& curve, Path& aPath,
                            int aIdx) const
            {
                typedef Point_2<Kernel>             Point;
                typedef CGAL::Qt::Converter<Kernel> Converter;
                Converter convert;

                Point ps(CGAL::to_double(curve[0].source().x()),
                         CGAL::to_double(curve[0].source().y()));
                Point pt(CGAL::to_double(curve[curve.number_of_subcurves()-1].target().x()),
                         CGAL::to_double(curve[curve.number_of_subcurves()-1].target().y()));

                if(aIdx == 0) aPath.moveTo(convert(ps));
                else aPath.lineTo(convert(ps));
                aPath.lineTo(convert(pt));

                //typedef Simple_cartesian<double> Linear_kernel;
                //commenting it gives errors
         /*       typedef Point_2<Kernel>             Linear_point;
                typedef CGAL::Qt::Converter<Kernel> Converter;
                Converter convert;

                Linear_point ps(CGAL::to_double(curve.source().x()),
                                CGAL::to_double(curve.source().y()));
                Linear_point pt(CGAL::to_double(curve.target().x()),
                                CGAL::to_double(curve.target().y()));

                if(aIdx == 0) aPath.moveTo(convert(ps));
                else aPath.lineTo(convert(ps));
                aPath.lineTo(convert(pt));*/

            }
        };

        struct Draw_polyline_curve {
            template <typename Polyline_segment_2, class Path>
            void operator()(Polyline_segment_2 const& curve, Path& aPath, int aIdx) const {
                //typedef Simple_cartesian<double> Linear_kernel;
                //commenting it gives errors

                typedef Point_2 <Kernel> Point;
                typedef Qt::Converter<Kernel> Converter;
                Converter convert;

                Point ps(CGAL::to_double(curve[0].source().x()),
                         CGAL::to_double(curve[0].source().y()));
                Point pt(CGAL::to_double(curve[curve.number_of_subcurves()-1].target().x()),
                         CGAL::to_double(curve[curve.number_of_subcurves()-1].target().y()));

                if (aIdx == 0) aPath.moveTo(convert(ps));
                else aPath.lineTo(convert(ps));
                aPath.lineTo(convert(pt));
            }
        };

        template <typename Polyline_boundary_pieces>
        class Polyline_boundary_pieces_graphics_item :
                public Boundary_pieces_graphics_item<Polyline_boundary_pieces,
                        Draw_polyline_curve, Polyline_bbox>
        {
            typedef Boundary_pieces_graphics_item<Polyline_boundary_pieces,
            Draw_polyline_curve, Polyline_bbox> Base;

        public:
            Polyline_boundary_pieces_graphics_item(Polyline_boundary_pieces* aPieces) :
                    Base(aPieces)
            {}
        };

        template<class Polyline_set,class Gps_traits>
        class Polyline_set_graphics_item :
                public Piecewise_set_graphics_item<Polyline_set, Gps_traits,
                        Draw_polyline_X_monotone_curve,
                        Polyline_X_monotone_bbox>
        {

            typedef Piecewise_set_graphics_item<Polyline_set, Gps_traits,
            Draw_polyline_X_monotone_curve,
            Polyline_X_monotone_bbox> Base;
        public:
            Polyline_set_graphics_item(Polyline_set* aSet, Gps_traits Polyline_traits) :
                    Base(aSet,Polyline_traits)
            {}
        };

    } // namespace Qt
} // namespace CGAL

