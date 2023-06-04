// Copyright (c) 2018-2020  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_POLYHEDRON_H
#define CGAL_DRAW_POLYHEDRON_H

#include <CGAL/license/Polyhedron.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/draw_face_graph.h>
#include <CGAL/Random.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/polyhedron_renderer.h>

#include <QVulkanWindow>
#include <QVulkanInstance>
#include <QVulkanFunctions>
#include <QLoggingCategory>


namespace CGAL
{
    template<class PolyhedronTraits_3,
        class PolyhedronItems_3,
        template < class T, class I, class A>
    class T_HDS,
        class Alloc>
    class VulkanWindow : public QVulkanWindow {
    public:
        QVulkanWindowRenderer* createRenderer()  {
            return new Polyhedron_renderer<PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>(this, m_poly);
        }
        void set_geometry(const CGAL_POLY_TYPE& poly) {
            m_poly = poly;
        }
    private:
        CGAL_POLY_TYPE m_poly;
    };

    template<class PolyhedronTraits_3,
        class PolyhedronItems_3,
        template < class T, class I, class A>
    class T_HDS,
        class Alloc>
    void draw(const CGAL_POLY_TYPE& apoly,
        const char* title = "Polyhedron Basic Viewer",
        bool nofill = false)
    {
#if defined(CGAL_TEST_SUITE)
        bool cgal_test_suite = true;
#else
        bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

        if (!cgal_test_suite)
        {
            //printf("number of vertices: %d\n", num_vertices(apoly));
            //CGAL_POLY_TYPE::Facet_const_iterator fi = apoly.faces_begin();
            //do {
              //  CGAL_POLY_TYPE::Halfedge_around_facet_const_circulator vi = apoly.halfedges_begin();;
            //} while (++fi != apoly.faces_end());
            //auto point_pmap = get(CGAL::vertex_point, apoly);

            //CGAL::Qt::init_ogl_context(4,3);
            int argc = 1;
            const char* argv[2] = { "polyhedron_viewer", nullptr };
            QApplication app(argc, const_cast<char**>(argv));
            QLoggingCategory::setFilterRules(QStringLiteral("qt.vulkan=true"));
            QVulkanInstance inst;
            inst.setLayers(QByteArrayList() << "VK_LAYER_KHRONOS_validation");
            if (!inst.create()) {
                qFatal("Failed to create a Vulkan instance: %d", inst.errorCode());
            }
            VulkanWindow<PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc> w;
            w.set_geometry(apoly);
            w.setVulkanInstance(&inst);
            w.showMaximized();
            QVulkanWindowRenderer* renderer = w.createRenderer();
            //renderer->geoData(apoly);
            //renderer->
            //renderer->startNextFrame();
            //for (auto v : vertices(apoly)) {
              //
            //}
            std::vector<float> vData{};
            w.requestUpdate();
            //SimpleFaceGraphViewerQt  mainwindow(app.activeWindow(), apoly, title, nofill);
            //mainwindow.show();
            app.exec();
        }
    }

#undef CGAL_POLY_TYPE

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_POLYHEDRON_H
