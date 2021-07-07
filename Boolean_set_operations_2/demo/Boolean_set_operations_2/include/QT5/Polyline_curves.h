// Copyright (c) 2009  GeometryFactory Sarl (France).
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
//
// Author(s) : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//			   Kabir Kedia <kabirkedia0111@gmail.com>
#ifndef CGAL_POLYLINE_CURVES_H
#define CGAL_POLYLINE_CURVES_H

#include <iostream>

#include <CGAL/Qt/Converter.h>
#include <CGAL/Arr_polyline_traits_2.h>

#include "Typedefs.h"
#include "QT5/Piecewise_set_graphics_item.h"
#include "QT5/Boundary_pieces_graphics_item.h"

namespace CGAL {
    namespace Qt {
        struct Polyline_X_monotone_bbox {
            template<typename X_monotone_polyline_segment_2>
            CGAL::Bbox_2 operator()(X_monotone_polyline_segment_2 const &aC) const { return aC.bbox(); }
        };

        struct Polyline_bbox {
            template<class Polyline_segment_2>
            CGAL::Bbox_2 operator()(Polyline_segment_2 const &aC) const {
                double x_min = to_double(aC.source().x());
                double x_max = to_double(aC.target().x());
                double y_min = to_double(aC.source().y());
                double y_max = to_double(aC.target().y());
                if (x_min > x_max) {
                    std::swap(x_min, x_max);
                    std::swap(y_min, y_max);
                }

                if (y_min > y_max) std::swap(y_min, y_max);

                return Bbox_2(x_min, y_min, x_max, y_max);
            }
        };

        struct Draw_polyline_X_monotone_curve {
            template <typename X_monotone_polyline_segment_2, typename Path>
            void operator()(X_monotone_polyline_segment_2 const& curve, Path& aPath,
                            int aIdx) const {

                typedef Simple_cartesian<double> Kernel;
                typedef Point_2 <Kernel> Linear_point;
                typedef CGAL::Qt::Converter<Kernel> Converter;
                Converter convert;
                Linear_point ps(CGAL::to_double(curve.source().x()),
                                CGAL::to_double(curve.source().y()));
                Linear_point pt(CGAL::to_double(curve.target().x()),
                                CGAL::to_double(curve.target().y()));

                if (aIdx == 0) aPath.moveTo(convert(ps));
                else aPath.lineTo(convert(ps));
                aPath.lineTo(convert(pt));
            }
        };

        struct Draw_polyline_curve {
            template <typename Linear_segment_2, class Path>
            void operator()(Linear_segment_2 const& curve, Path& aPath, int aIdx) const {
                //commenting it gives errors
                typedef Simple_cartesian<double> Kernel;
                typedef Point_2 <Kernel> Linear_point;
                typedef Qt::Converter<Kernel> Converter;
                Converter convert;

                Linear_point ps(CGAL::to_double(curve.source().x()),
                                CGAL::to_double(curve.source().y()));
                Linear_point pt(CGAL::to_double(curve.target().x()),
                                CGAL::to_double(curve.target().y()));

                if (aIdx == 0) aPath.moveTo(convert(ps));
                else aPath.lineTo(convert(ps));
                aPath.lineTo(convert(pt));
            }
        };

        template<typename Polyline_boundary_pieces>
        class Polyline_boundary_pieces_graphics_item :
                public Boundary_pieces_graphics_item<Polyline_boundary_pieces,
                                                     Draw_polyline_curve, Polyline_bbox>
                                                     {
            typedef Boundary_pieces_graphics_item <Polyline_boundary_pieces,
                                                   Draw_polyline_curve,
                                                   Polyline_bbox> Base;
        public :
            Polyline_boundary_pieces_graphics_item(Polyline_boundary_pieces *aPieces) :
                    Base(aPieces)
                    {}
        };

        template<typename Polyline_set, typename Gps_traits>
        class Polyline_set_graphics_item :
                public Piecewise_set_graphics_item<Polyline_set, Gps_traits,
                                                   Draw_polyline_X_monotone_curve,
                                                   Polyline_X_monotone_bbox> {

            typedef Piecewise_set_graphics_item <Polyline_set, Gps_traits,
                                                 Draw_polyline_X_monotone_curve,
                                                 Polyline_X_monotone_bbox> Base;

        public:

            Polyline_set_graphics_item(Polyline_set *aSet, Gps_traits Polyline_traits) :
                    Base(aSet, Polyline_traits) {}
        };
    }//namespace Qt
}//namespace CGAL

#endif //CGAL_POLYLINE_CURVES_H
